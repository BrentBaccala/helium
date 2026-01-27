#!/bin/bash

# Worker Log Analysis Script
# Author: qwen3-coder
# Usage: ./analyze_workers.sh [log_file]

LOG_FILE="${1:-screenlog.0}"

if [[ ! -f "$LOG_FILE" ]]; then
    echo "Error: Log file '$LOG_FILE' not found."
    exit 1
fi

echo "Analyzing worker processes in: $LOG_FILE"
echo "=========================================="

# Extract all data once to avoid repeated grep calls
mapfile -t START_LOGS < <(grep '^Started worker' "$LOG_FILE")
mapfile -t EXIT_LOGS < <(grep '^Worker [0-9]' "$LOG_FILE" | grep -E 'exiting|exited')

# Parse started PIDs
declare -A started_pids
for line in "${START_LOGS[@]}"; do
    pid=$(echo "$line" | awk '{gsub(/\./,"",$3); print $3}')
    started_pids["$pid"]="$line"
done

# Parse exited PIDs
declare -A exited_pids
for line in "${EXIT_LOGS[@]}"; do
    pid=$(echo "$line" | awk '{print $2}')
    exited_pids["$pid"]="$line"
done

# Convert associative array keys to indexed arrays for processing
all_started=("${!started_pids[@]}")
all_exited=("${!exited_pids[@]}")

# Find PIDs that started but never exited (may still be running)
non_exited_pids=()
for pid in "${all_started[@]}"; do
    if [[ -z "${exited_pids[$pid]+isset}" ]]; then
        non_exited_pids+=("$pid")
    fi
done

# Find PIDs that exited but never had a start message (shouldn't happen in clean logs)
orphaned_exit_pids=()
for pid in "${all_exited[@]}"; do
    if [[ -z "${started_pids[$pid]+isset}" ]]; then
        orphaned_exit_pids+=("$pid")
    fi
done

echo "SUMMARY:"
echo "--------"
echo "Total started: ${#all_started[@]}"
echo "Total exited: ${#all_exited[@]}"
echo "Started but not yet logged as exited: ${#non_exited_pids[@]}"
if [[ ${#orphaned_exit_pids[@]} -gt 0 ]]; then
    echo "ERROR: PIDs with exit logs but no start logs: ${#orphaned_exit_pids[@]}"
fi
echo

echo "DETAILED ANALYSIS:"
echo "------------------"

# Categorize non-exited PIDs
running_pids=()
zombie_pids=()
vanished_pids=()
for pid in "${non_exited_pids[@]}"; do
    # Check if process exists and its state
    if ps -p "$pid" -o stat= 2>/dev/null 1>&2; then
        state=$(ps -p "$pid" -o stat= 2>/dev/null)
        # Remove leading/trailing whitespace
        state=$(echo "$state" | xargs)
        if [[ "$state" == *"Z"* ]]; then
            zombie_pids+=("$pid")
        else
            running_pids+=("$pid")
        fi
    else
        # Process doesn't exist
        vanished_pids+=("$pid")
    fi
done

# Display running processes
if [[ ${#running_pids[@]} -gt 0 ]]; then
    echo "RUNNING PROCESSES (started but not yet exited):"
    for pid in "${running_pids[@]}"; do
        echo "  PID $pid: ${started_pids[$pid]}"
    done
    echo
else
    echo "RUNNING PROCESSES: None found"
    echo
fi

# Display zombie processes
if [[ ${#zombie_pids[@]} -gt 0 ]]; then
    echo "ZOMBIE PROCESSES (started but not yet exited):"
    for pid in "${zombie_pids[@]}"; do
        echo "  PID $pid: ${started_pids[$pid]}"
    done
    echo
else
    echo "ZOMBIE PROCESSES: None found"
    echo
fi

# Display vanished processes
if [[ ${#vanished_pids[@]} -gt 0 ]]; then
    echo "VANISHED PROCESSES (started but not running, no exit log):"
    for pid in "${vanished_pids[@]}"; do
        echo "  PID $pid: ${started_pids[$pid]}"
    done
    echo
else
    echo "VANISHED PROCESSES: None found"
    echo
fi

if [[ ${#orphaned_exit_pids[@]} -gt 0 ]]; then
    echo "ORPHANED EXITS (no matching start log):"
    for pid in "${orphaned_exit_pids[@]}"; do
        echo "  PID $pid: ${exited_pids[$pid]}"
    done
    echo
fi

echo "FINAL COUNTS:"
echo "-------------"
echo "Normal (running): ${#running_pids[@]}"
echo "Zombies: ${#zombie_pids[@]}"
echo "Completed (logged exit): ${#all_exited[@]}"
if [[ ${#vanished_pids[@]} -gt 0 ]]; then
    echo "Vanished (started but missing): ${#vanished_pids[@]}"
fi
if [[ ${#orphaned_exit_pids[@]} -gt 0 ]]; then
    echo "Orphaned (exit log without start): ${#orphaned_exit_pids[@]}"
fi
