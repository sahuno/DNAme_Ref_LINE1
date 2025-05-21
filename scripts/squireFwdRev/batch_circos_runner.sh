#!/usr/bin/env bash
# batch_circos_runner.sh
# Usage: ./batch_circos_runner.sh --input-dir <tsv_directory> --outdir <output_directory> --script <path_to_circos_plot.R>
# Usage: bash /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/scripts/squireFwdRev/batch_circos_runner.sh --input-dir \
# /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/outputs/LocusSpecRepeatDE_lowQuaSamplesDropped/upsetPlots/separateHistoneModifications/QSTATi/nonSmallRNA \
# --outdir /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/outputs/LocusSpecRepeatDE_lowQuaSamplesDropped/circosPlots \
# --script /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/scripts/squireFwdRev/plotCircos.R --shared FALSE


# ./batch_circos_runner.sh --input-dir <tsv_directory> --outdir <output_directory> --script <path_to_circos_plot.R>
# Usage: bash /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/scripts/squireFwdRev/batch_circos_runner.sh \
# --input-dir /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/outputs/LocusSpecRepeatDE_lowQuaSamplesDropped/upsetPlots/separateHistoneModifications/SETDB1i/nonSmallRNA \
# --outdir /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/outputs/LocusSpecRepeatDE_lowQuaSamplesDropped/circosPlots/SETDB1i/ \
# --script /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/scripts/squireFwdRev/plotCircos.R --shared FALSE


####shared repeats; QSTATi
# bash /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/scripts/squireFwdRev/batch_circos_runner.sh \
# --input-dir /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/outputs/LocusSpecRepeatDE_lowQuaSamplesDropped/upsetPlots/separateHistoneModifications/QSTATi/nonSmallRNA \
# --outdir /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/outputs/LocusSpecRepeatDE_lowQuaSamplesDropped/circosPlots/QSTATi/sharedRepeats/ \
# --script /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/scripts/squireFwdRev/plotCircos.R --shared

####shared repeats: SETDB1i
# ./batch_circos_runner.sh --input-dir <tsv_directory> --outdir <output_directory> --script <path_to_circos_plot.R>
# bash /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/scripts/squireFwdRev/batch_circos_runner.sh \
# --input-dir /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/outputs/LocusSpecRepeatDE_lowQuaSamplesDropped/upsetPlots/separateHistoneModifications/SETDB1i/nonSmallRNA \
# --outdir /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/outputs/LocusSpecRepeatDE_lowQuaSamplesDropped/circosPlots/SETDB1i/sharedRepeats/ \
# --script /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/scripts/squireFwdRev/plotCircos.R --shared





# /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/outputs/LocusSpecRepeatDE_lowQuaSamplesDropped/upsetPlots/separateHistoneModifications/SETDB1i/nonSmallRNA
# batch_circos_runner.sh
# Usage: $0 --input-dir <directory> --outdir <directory> --script <circos_plot.R> [--shared]
set -euo pipefail

print_usage() {
  echo "Usage: $0 --input-dir <directory> --outdir <directory> --script <circos_plot.R> [--shared]"
  exit 1
}

# defaults
SHARED_FLAG=false

# parse args
while [[ $# -gt 0 ]]; do
  case $1 in
    --input-dir)
      INPUT_DIR="$2"; shift 2;;
    --outdir)
      OUTDIR="$2"; shift 2;;
    --script)
      PLOT_SCRIPT="$2"; shift 2;;
    --shared)
      SHARED_FLAG=true; shift;;
    *)
      echo "Unknown argument: $1" >&2
      print_usage;;
  esac
done

# validate
if [[ -z "${INPUT_DIR:-}" || -z "${OUTDIR:-}" || -z "${PLOT_SCRIPT:-}" ]]; then
  echo "Error: --input-dir, --outdir, and --script are required." >&2
  print_usage
fi
if [[ ! -d "$INPUT_DIR" ]]; then
  echo "Error: input-dir '$INPUT_DIR' does not exist." >&2
  exit 2
fi

# create outdir
mkdir -p "$OUTDIR"

# loop over TSV files
for tsv in "$INPUT_DIR"/*.tsv; do
  fname=$(basename "$tsv" .tsv)
  keyword="$fname"
  echo "Processing: $tsv  (keyword = $keyword, shared = $SHARED_FLAG)"

  # build Rscript arguments
  args=(--input "$tsv" --outdir "$OUTDIR" --keyword "$keyword")
  if [[ "$SHARED_FLAG" == "true" ]]; then
    args+=(--shared)
  fi

  # try plotting; on error just warn & continue
  if ! Rscript "$PLOT_SCRIPT" "${args[@]}"; then
    echo "Warning: plotting failed for '$fname' (likely no data after --shared filter); skipping."
    continue
  fi
done

echo "Done.  Generated plots (skipping failures) in: $OUTDIR"
