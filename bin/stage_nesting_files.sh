#!/bin/bash
# Stage offline nesting archive files for a given run period.
# Copies files directly if available, otherwise extracts from .tar.gz archives.
# Stages DATE_START-1 through DATE_END+1 so HYCOM can interpolate across the
# full run period including the first and last time steps.
#
# Usage: stage_nesting_files.sh DIR_NST WORK_NST DATE_START DATE_END [OPTIONS]
#
#   DIR_NST          Source directory containing nesting archive files or tar_files/
#   WORK_NST         Destination directory (created if it does not exist)
#   DATE_START       First day of run period (YYYY-MM-DD)
#   DATE_END         Last day of run period (YYYY-MM-DD)
#   --no-fabm        Skip BGC (FABM) archive files
#   --skip-existing  Skip dates where both .a and .b files already exist in WORK_NST

usage() {
    echo "Usage: $0 DIR_NST WORK_NST DATE_START DATE_END [--no-fabm] [--skip-existing]"
    exit 1
}

[ $# -lt 4 ] && usage

DIR_NST=$1
WORK_NST=$2
DATE_START=$3
DATE_END=$4
WITH_FABM=true
SKIP_EXISTING=false
for arg in "${@:5}"; do
    case "$arg" in
        --no-fabm)       WITH_FABM=false ;;
        --skip-existing) SKIP_EXISTING=true ;;
        *) echo "Unknown option: $arg"; usage ;;
    esac
done

mkdir -p "${WORK_NST}" || { echo "Cannot create ${WORK_NST}"; exit 1; }
cd "${WORK_NST}"       || { echo "Cannot cd to ${WORK_NST}";   exit 1; }

d_start=$(date -u -d "$DATE_START - 1 day" +%s)
d_end=$(date -u -d "$DATE_END + 1 day" +%s)

# Build the set of needed filenames (for --skip-existing and direct-copy paths)
declare -A needed_physics needed_fabm
for (( d=d_start; d<=d_end; d+=86400 )); do
    DATE_NOW=$(date -u -d "@$d" +%Y-%m-%d)
    YYYY=$(date -d "$DATE_NOW" +%Y)
    DOY=$(date -d  "$DATE_NOW" +%j)
    afile=archv.${YYYY}_${DOY}_00.a
    bfile=archv.${YYYY}_${DOY}_00.b
    afile_fabm=archv_fabm.${YYYY}_${DOY}_00.a
    bfile_fabm=archv_fabm.${YYYY}_${DOY}_00.b

    # Skip if already present
    if $SKIP_EXISTING && [ -f "$afile" ] && [ -f "$bfile" ]; then
        if ! $WITH_FABM || ( [ -f "$afile_fabm" ] && [ -f "$bfile_fabm" ] ); then continue; fi
    fi

    # Copy directly if individual files are available
    if [ -f "${DIR_NST}/${afile}" ]; then
        cp "${DIR_NST}/${afile}" "${DIR_NST}/${bfile}" .
        if $WITH_FABM && [ -f "${DIR_NST}/${afile_fabm}" ]; then
            cp "${DIR_NST}/${afile_fabm}" "${DIR_NST}/${bfile_fabm}" .
        fi
        continue
    fi

    # Otherwise queue for tar extraction
    needed_physics["$afile"]=1
    needed_physics["$bfile"]=1
    if $WITH_FABM; then
        needed_fabm["$afile_fabm"]=1
        needed_fabm["$bfile_fabm"]=1
    fi
done

# Extract all needed files from tar archives — one pass per archive
if [ -d "${DIR_NST}/tar_files" ] && [ ${#needed_physics[@]} -gt 0 ]; then
    for f in "${DIR_NST}/tar_files/archv."*.tar.gz; do
        [ -f "$f" ] || continue
        # Collect files from this archive that we actually need
        to_extract=()
        while IFS= read -r member; do
            [[ -v needed_physics["$member"] ]] && to_extract+=("$member")
        done < <(tar -tzf "$f" 2>/dev/null)
        if [ ${#to_extract[@]} -gt 0 ]; then
            echo "Extracting ${#to_extract[@]} files from $(basename "$f")"
            tar -xzf "$f" "${to_extract[@]}"
        fi
    done
fi

if $WITH_FABM && [ -d "${DIR_NST}/tar_files" ] && [ ${#needed_fabm[@]} -gt 0 ]; then
    for f in "${DIR_NST}/tar_files/archv_fabm."*.tar.gz; do
        [ -f "$f" ] || continue
        to_extract=()
        while IFS= read -r member; do
            [[ -v needed_fabm["$member"] ]] && to_extract+=("$member")
        done < <(tar -tzf "$f" 2>/dev/null)
        if [ ${#to_extract[@]} -gt 0 ]; then
            echo "Extracting ${#to_extract[@]} files from $(basename "$f")"
            tar -xzf "$f" "${to_extract[@]}"
        fi
    done
fi
