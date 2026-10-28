"""Read a tabix-indexed, bgzip-compressed summary-statistics file over HTTP range requests.

Why this exists: the FinnGen R12 release publishes ~815 MB per endpoint plus a ~1.7 MB
`.tbi` index. Downloading eight of them for a 600-variant Mendelian-randomisation analysis
costs 6.5 GB; the index plus ~600 single-position queries costs about 30 MB per endpoint.
Throughput measured on this project's connection is 0.05-0.1 MB/s per socket, which makes
the difference between minutes and hours.

The GWAS files this project downloads in full are read with pandas, not with this reader; it
is deliberately narrow: bgzip blocks, the .tbi binning/linear index, and generic (tabix -p
vcf-free) text files. No third-party dependency, because neither pysam nor the tabix binary
is installed here and htslib's remote support is not guaranteed on Windows.

Usage:
    from _tabix_http import TabixHttp
    t = TabixHttp("https://.../finngen_R12_I9_AF.gz")
    for line in t.query("1", 154426264, 154426264):
        ...
    t.header                      # the '#chrom pos ref ...' line, as a list of column names
"""
import concurrent.futures as cf
import gzip
import struct
import threading
import time
import urllib.request
import zlib

UA = {"User-Agent": "Mozilla/5.0 (research; tabix range reader)"}


def _reg2bins(beg, end):
    """htslib reg2bins: the index bins a 0-based half-open interval can fall in."""
    if beg >= end:
        return []
    if end >= 1 << 29:
        end = 1 << 29
    end -= 1
    bins = [0]
    bins += range(1 + (beg >> 26), 2 + (end >> 26))
    bins += range(9 + (beg >> 23), 10 + (end >> 23))
    bins += range(73 + (beg >> 20), 74 + (end >> 20))
    bins += range(585 + (beg >> 17), 586 + (end >> 17))
    bins += range(4681 + (beg >> 14), 4682 + (end >> 14))
    return bins


class TabixHttp:
    """Minimal tabix client for one remote bgzip file + its local .tbi index."""

    def __init__(self, url, tbi_path, timeout=120, retries=5):
        self.url = url
        self.timeout = timeout
        self.retries = retries
        self._blocks = {}
        self._lock = threading.Lock()
        self._load_index(str(tbi_path))
        self.header = self._read_header()

    # ---------------------------------------------------------------- transport
    def _get(self, start, end):
        """Bytes [start, end] inclusive, with retries. Servers here answer HEAD with 404."""
        for attempt in range(self.retries):
            try:
                req = urllib.request.Request(
                    self.url, headers={**UA, "Range": f"bytes={start}-{end}"})
                with urllib.request.urlopen(req, timeout=self.timeout) as r:
                    if r.status not in (200, 206):
                        raise RuntimeError(f"status {r.status}")
                    return r.read()
            except Exception:
                if attempt == self.retries - 1:
                    raise
                time.sleep(1.0 + 1.5 * attempt)
        raise RuntimeError("unreachable")

    # ---------------------------------------------------------------- bgzf
    def block(self, coffset):
        """Decompress the bgzip block starting at `coffset`; return (plain bytes, next offset).

        The whole block is pulled with one request: the BSIZE field lives at byte 16, but
        asking for the largest possible block (64 KB 鈥?BSIZE is a uint16, so a bgzip block
        can never exceed it) and trimming costs one round trip instead of two, and on a
        high-latency link the round trips dominate.
        """
        with self._lock:
            hit = self._blocks.get(coffset)
        if hit is not None:
            return hit
        buf = self._get(coffset, coffset + (1 << 16) - 1)
        if len(buf) < 18 or buf[0] != 31 or buf[1] != 139:
            raise ValueError(f"no bgzip block at {coffset}")
        xlen = struct.unpack("<H", buf[10:12])[0]
        extra = buf[12:12 + xlen]
        bsize = None
        i = 0
        while i + 4 <= len(extra):
            si1, si2, slen = extra[i], extra[i + 1], struct.unpack("<H", extra[i + 2:i + 4])[0]
            if si1 == 66 and si2 == 67:                       # 'BC' subfield carries BSIZE
                bsize = struct.unpack("<H", extra[i + 4:i + 6])[0]
            i += 4 + slen
        if bsize is None:
            raise ValueError(f"block at {coffset} has no BSIZE")
        total = bsize + 1
        if len(buf) < total:
            buf = self._get(coffset, coffset + total - 1)
        if len(buf) < total:
            raise ValueError(f"short read at {coffset}: {len(buf)} < {total}")
        payload = buf[12 + xlen:total - 8]
        plain = zlib.decompress(payload, -15)
        out = (plain, coffset + total)
        with self._lock:
            self._blocks[coffset] = out
        return out

    def _slice(self, vbeg, vend):
        """Plain text covering [vbeg, vend), decompressing only the blocks in between."""
        start_coff, end_coff = vbeg >> 16, vend >> 16
        coff, uoff = start_coff, vbeg & 0xFFFF
        parts = []
        while coff <= end_coff:
            plain, nxt = self.block(coff)
            lo = uoff if coff == start_coff else 0
            hi = (vend & 0xFFFF) if coff == end_coff else len(plain)
            hi = min(hi, len(plain))
            if lo < hi:
                parts.append(plain[lo:hi])
            if nxt <= coff:                                   # bgzip EOF marker
                break
            coff, uoff = nxt, 0
        return b"".join(parts)

    # ---------------------------------------------------------------- index
    def _load_index(self, tbi_path):
        with open(tbi_path, "rb") as f:
            raw = f.read()
        data = gzip.decompress(raw) if raw[:2] == b"\x1f\x8b" else raw
        if data[:4] != b"TBI\x01":
            raise ValueError("not a tabix index")
        (n_ref, self.fmt, self.col_seq, self.col_beg, self.col_end,
         self.meta, self.skip, l_nm) = struct.unpack("<8i", data[4:36])
        p = 36
        names = data[p:p + l_nm].split(b"\x00")
        p += l_nm
        self.refs = {}
        for name in names[:n_ref]:
            n_bin = struct.unpack("<i", data[p:p + 4])[0]
            p += 4
            bins = {}
            for _ in range(n_bin):
                bin_id, n_chunk = struct.unpack("<Ii", data[p:p + 8])
                p += 8
                chunks = [struct.unpack("<QQ", data[p + 16 * k:p + 16 * k + 16])
                          for k in range(n_chunk)]
                p += 16 * n_chunk
                bins[bin_id] = chunks
            n_intv = struct.unpack("<i", data[p:p + 4])[0]
            p += 4
            linear = [struct.unpack("<Q", data[p + 8 * k:p + 8 * k + 8])[0]
                      for k in range(n_intv)]
            p += 8 * n_intv
            self.refs[name.decode()] = (bins, linear)

    def _read_header(self):
        plain, _ = self.block(0)
        first = plain.split(b"\n")[0].decode()
        return first.lstrip("#").split("\t")

    # ---------------------------------------------------------------- query
    def query(self, chrom, start, end):
        """Lines overlapping the 1-based inclusive interval, as lists of fields."""
        c = str(chrom)
        if c.endswith(".0"):                    # chromosome columns read as floats
            c = c[:-2]
        ref = self.refs.get(c)
        if ref is None:
            ref = self.refs.get(c[3:] if c.lower().startswith("chr") else "chr" + c)
        if ref is None:
            return []
        bins, linear = ref
        beg, stop = start - 1, end                          # 0-based half-open
        min_off = linear[beg >> 14] if (beg >> 14) < len(linear) else 0
        chunks = []
        for b in _reg2bins(beg, stop):
            for c in bins.get(b, ()):
                if c[1] > min_off:
                    chunks.append(c)
        if not chunks:
            return []
        chunks.sort()
        merged = [list(chunks[0])]
        for c in chunks[1:]:
            if c[0] <= merged[-1][1]:
                merged[-1][1] = max(merged[-1][1], c[1])
            else:
                merged.append(list(c))

        out = []
        for vbeg, vend in merged:
            for line in self._slice(vbeg, vend).split(b"\n"):
                if not line:
                    continue
                s = line.decode()
                if s.startswith("#"):
                    continue
                f = s.split("\t")
                try:
                    pos = int(f[self.col_beg - 1])
                except (IndexError, ValueError):
                    continue
                if start <= pos <= end:
                    out.append(f)
        return out


def fetch_worker(reader, jobs):
    rows = []
    for snp, chrom, pos in jobs:
        for f in reader.query(str(chrom), int(pos), int(pos)):
            rows.append((snp, pos, f))
    return rows


def fetch(reader, targets, threads=8):
    """Run `targets` ([(snp, chrom, pos)]) against the reader on a small thread pool."""
    batches = [targets[i:i + 25] for i in range(0, len(targets), 25)]
    rows = []
    with cf.ThreadPoolExecutor(max_workers=threads) as ex:
        for r in ex.map(lambda b: fetch_worker(reader, b), batches):
            rows.extend(r)
    return rows

