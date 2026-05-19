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

for (( d=$(date -u -d "$DATE_START - 1 day" +%s); d<=$(date -u -d "$DATE_END + 1 day" +%s); d+=86400 )); do
    DATE_NOW=$(date -u -d "@$d" +%Y-%m-%d)
    YYYY=$(date -d "$DATE_NOW" +%Y)
    DOY=$(date -d  "$DATE_NOW" +%j)

    afile=archv.${YYYY}_${DOY}_00.a
    bfile=archv.${YYYY}_${DOY}_00.b
    afile_fabm=archv_fabm.${YYYY}_${DOY}_00.a
    bfile_fabm=archv_fabm.${YYYY}_${DOY}_00.b

    if $SKIP_EXISTING && [ -f "$afile" ] && [ -f "$bfile" ]; then continue; fi

    # copy directly if files are available individually
    if [ -f "${DIR_NST}/${afile}" ]; then
        cp "${DIR_NST}/${afile}" "${DIR_NST}/${bfile}" .
        if $WITH_FABM && [ -f "${DIR_NST}/${afile_fabm}" ]; then
            cp "${DIR_NST}/${afile_fabm}" "${DIR_NST}/${bfile_fabm}" .
        fi
        continue
    fi

    # otherwise extract from tar archives (files may be grouped by DOY range)
    if [ -d "${DIR_NST}/tar_files" ]; then
        for f in "${DIR_NST}/tar_files/archv.${YYYY}_"*.tar.gz; do
            range=$(echo "${f#${DIR_NST}/tar_files/archv.${YYYY}_}" | sed 's/\.tar\.gz//')
            start=$((10#${range%_*}))
            end=$((10#${range#*_}))
            if (( 10#$DOY >= start && 10#$DOY <= end )); then
                tar -xzf "${DIR_NST}/tar_files/archv.${YYYY}_${range}.tar.gz" "$afile" "$bfile"
            fi
        done
        if $WITH_FABM; then
            for f in "${DIR_NST}/tar_files/archv_fabm.${YYYY}_"*.tar.gz; do
                range=$(echo "${f#${DIR_NST}/tar_files/archv_fabm.${YYYY}_}" | sed 's/\.tar\.gz//')
                start=$((10#${range%_*}))
                end=$((10#${range#*_}))
                if (( 10#$DOY >= start && 10#$DOY <= end )); then
                    tar -xzf "${DIR_NST}/tar_files/archv_fabm.${YYYY}_${range}.tar.gz" "$afile_fabm" "$bfile_fabm"
                fi
            done
        fi
    fi
done
