import argparse
from importlib.resources import files

from eber_scope import __version__
from eber_scope.detect import CHEMISTRIES, run_detection
from eber_scope.merge_calls import merge_calls


def build_parser():
    parser = argparse.ArgumentParser(prog="eber-scope")
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    subparsers = parser.add_subparsers(dest="command", required=True)

    detect = subparsers.add_parser("detect", help="Detect EBER1-positive cell barcodes")
    detect.add_argument("--r1", required=True)
    detect.add_argument("--r2", required=True)
    detect.add_argument("--sample", required=True)
    detect.add_argument("--output-prefix", required=True)
    detect.add_argument("--reference", default=None)
    detect.add_argument("--barcode-whitelist")
    detect.add_argument("--chemistry", choices=sorted(CHEMISTRIES), default="5p-v3")
    detect.add_argument("--left-anchor", default="CTTCCGATCT")
    detect.add_argument("--right-anchor", default="TTTCTTATATGGG")
    detect.add_argument("--anchor-mismatches", type=int, default=1)
    detect.add_argument("--threshold", type=float, default=0.60)
    detect.add_argument("--max-reads", type=int)

    merge = subparsers.add_parser("merge", help="Combine targeted and existing cell calls")
    merge.add_argument("--existing", required=True)
    merge.add_argument("--targeted", required=True)
    merge.add_argument("--output", required=True)
    return parser


def main(argv=None):
    args = build_parser().parse_args(argv)
    if args.command == "detect":
        reference = args.reference
        if reference is None:
            reference = files("eber_scope").joinpath("references/EBER1_NC_007605.1.fa")
        run_detection(
            r1_path=args.r1,
            r2_path=args.r2,
            output_prefix=args.output_prefix,
            sample=args.sample,
            reference_path=reference,
            barcode_whitelist=args.barcode_whitelist,
            chemistry=args.chemistry,
            left_anchor=args.left_anchor,
            right_anchor=args.right_anchor,
            anchor_mismatches=args.anchor_mismatches,
            threshold=args.threshold,
            max_reads=args.max_reads,
        )
    else:
        merge_calls(args.existing, args.targeted, args.output)
