#!/usr/bin/env python3
import argparse
import gzip

def open_maybe_gz(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")

def fasta_records(path):
    """Yield (header_without_>, sequence_string)"""
    with open_maybe_gz(path) as f:
        header = None
        seq_chunks = []
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_chunks)
                header = line[1:]
                seq_chunks = []
            else:
                seq_chunks.append(line)
        if header is not None:
            yield header, "".join(seq_chunks)

def wrap_fasta(seq, width=60):
    for i in range(0, len(seq), width):
        yield seq[i:i+width]

def main():
    ap = argparse.ArgumentParser(
        description="Extract region from FASTA (header-match mode or contig slice mode)."
    )
    ap.add_argument("-i", "--input", required=True, help="Input FASTA (.fa/.fa.gz)")
    ap.add_argument("-o", "--output", required=True, help="Output FASTA")
    ap.add_argument("--key", default=None, help="Mode 1: exact header after '>' (e.g. 'chrR:9307-45306')")
    ap.add_argument("-c", "--contig", default=None, help="Mode 2: contig name (e.g. 'chrR')")
    ap.add_argument("-s", "--start", type=int, default=None, help="Mode 2: start (1-based, inclusive)")
    ap.add_argument("-e", "--end", type=int, default=None, help="Mode 2: end (1-based, inclusive)")
    args = ap.parse_args()

    out_lines = []

    if args.key is not None:
        # Mode 1: exact header match
        found = False
        for header, seq in fasta_records(args.input):
            if header == args.key:
                out_lines.append(">" + header)
                out_lines.extend(wrap_fasta(seq))
                found = True
                break
        if not found:
            raise SystemExit(f"ERROR: header '{args.key}' not found in FASTA.")
        with open(args.output, "w") as w:
            w.write("\n".join(out_lines) + "\n")
        print(f"Done: extracted record with header >{args.key} -> {args.output}")
        return

    # Mode 2: slice from contig by coordinates
    if args.contig is None or args.start is None or args.end is None:
        raise SystemExit("ERROR: Provide either --key (mode 1) or --contig/--start/--end (mode 2).")

    contig = args.contig
    start = args.start
    end = args.end

    if start < 1 or end < 1 or end < start:
        raise SystemExit("ERROR: Require 1-based coordinates with end >= start.")

    # assume coordinates are on concatenated contig sequence (1-based)
    seq = None
    for header, s in fasta_records(args.input):
        # match "header starts with contig" (common: 'chrR' or 'chrR something')
        # but also accept exact header == contig
        first_field = header.split()[0]
        if first_field == contig:
            seq = s
            matched_header = header
            break

    if seq is None:
        raise SystemExit(f"ERROR: contig '{contig}' not found in FASTA (matched by first field).")

    if end > len(seq):
        end = len(seq)

    subseq = seq[start-1:end]  # 1-based -> 0-based slice
    out_header = f"{contig}:{start}-{end}"

    out_lines.append(">" + out_header)
    out_lines.extend(wrap_fasta(subseq))
    with open(args.output, "w") as w:
        w.write("\n".join(out_lines) + "\n")

    print(f"Done: extracted {out_header} from '{matched_header}' -> {args.output}")

if __name__ == "__main__":
    main()
