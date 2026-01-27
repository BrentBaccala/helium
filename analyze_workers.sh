#!/bin/bash

# Worker Log Analysis Script
# Author: qwen3-coder
# Usage: ./analyze_workers.sh [log_file] [--debug]

# Parse command line arguments
DEBUG=false
LOG_FILE="screenlog.0"

for arg in "$@"; do
    case $arg in
        --debug)
            DEBUG=true
            ;;
        *)
            # If it's not --debug and doesn't start with -, assume it's the log file
            if [[ $arg != -* ]]; then
                LOG_FILE="$arg"
            fi
            ;;
    esac
done

if [[ ! -f "$LOG_FILE" ]]; then
    echo "Error: Log file '$LOG_FILE' not found."
    exit 1
fi

if [[ "$DEBUG" == true ]]; then
    echo "[DEBUG] Using log file: $LOG_FILE"
    echo "[DEBUG] Debug mode enabled"
    echo
fi

echo "Analyzing worker processes in: $LOG_FILE"
echo "=========================================="

# Extract all data once to avoid repeated grep calls
if [[ "$DEBUG" == true ]]; then
    echo "[DEBUG] Running: grep -a '^Started worker' '$LOG_FILE'"
fi
mapfile -t START_LOGS < <(grep -a '^Started worker' "$LOG_FILE")
if [[ "$DEBUG" == true ]]; then
    echo "[DEBUG] Found ${#START_LOGS[@]} start logs:"
    for i in "${!START_LOGS[@]}"; do
        echo "  [DEBUG] Start log $i: ${START_LOGS[$i]}"
    done
    echo
fi

if [[ "$DEBUG" == true ]]; then
    echo "[DEBUG] Running: grep -a '^Worker [0-9]' '$LOG_FILE' | grep -a -E 'exiting|exited'"
fi
mapfile -t EXIT_LOGS < <(grep -a '^Worker [0-9]' "$LOG_FILE" | grep -a -E 'exiting|exited')
if [[ "$DEBUG" == true ]]; then
    echo "[DEBUG] Found ${#EXIT_LOGS[@]} exit logs:"
    for i in "${!EXIT_LOGS[@]}"; do
        echo "  [DEBUG] Exit log $i: ${EXIT_LOGS[$i]}"
    done
    echo
fi

# Parse started PIDs
declare -A started_pids
if [[ "$DEBUG" == true ]]; then
    echo "[DEBUG] Parsing started PIDs..."
fi
for line in "${START_LOGS[@]}"; do
    if [[ "$DEBUG" == true ]]; then
        echo "  [DEBUG] Processing start line: $line"
    fi
    pid=$(echo "$line" | awk '{gsub(/\./,"",$3); print $3}')
    if [[ "$DEBUG" == true ]]; then
        echo "    [DEBUG] Extracted PID: $pid"
    fi
    started_pids["$pid"]="$line"
done

if [[ "$DEBUG" == true ]]; then
    echo "[DEBUG] Started PIDs array contents:"
    for pid in "${!started_pids[@]}"; do
        echo "  [DEBUG] PID $pid -> ${started_pids[$pid]}"
    done
    echo
fi

# Parse exited PIDs
declare -A exited_pids
if [[ "$DEBUG" == true ]]; then
    echo "[DEBUG] Parsing exited PIDs..."
fi
for line in "${EXIT_LOGS[@]}"; do
    if [[ "$DEBUG" == true ]]; then
        echo "  [DEBUG] Processing exit line: $line"
    fi
    pid=$(echo "$line" | awk '{print $2}')
    if [[ "$DEBUG" == true ]]; then
        echo "    [DEBUG] Extracted PID: $pid"
    fi
    exited_pids["$pid"]="$line"
done

if [[ "$DEBUG" == true ]]; then
    echo "[DEBUG] Exited PIDs array contents:"
    for pid in "${!exited_pids[@]}"; do
        echo "  [DEBUG] PID $pid -> ${exited_pids[$pid]}"
    done
    echo
fi

# Convert associative array keys to indexed arrays for processing
all_started=("${!started_pids[@]}")
all_exited=("${!exited_pids[@]}")

if [[ "$DEBUG" == true ]]; then
    echo "[DEBUG] All started PIDs: ${all_started[*]}"
    echo "[DEBUG] All exited PIDs: ${all_exited[*]}"
    echo
fi

# Find PIDs that started but never exited (may still be running)
non_exited_pids=()
if [[ "$DEBUG" == true ]]; then
    echo "[DEBUG] Finding PIDs that started but never exited..."
fi
for pid in "${all_started[@]}"; do
    if [[ "$DEBUG" == true ]]; then
        echo "  [DEBUG] Checking if PID $pid has an exit log..."
    fi
    if [[ -z "${exited_pids[$pid]+isset}" ]]; then
        if [[ "$DEBUG" == true ]]; then
            echo "    [DEBUG] PID $pid has NO exit log, adding to non_exited_pids"
        fi
        non_exited_pids+=("$pid")
    else
        if [[ "$DEBUG" == true ]]; then
            echo "    [DEBUG] PID $pid HAS an exit log"
        fi
    fi
done

if [[ "$DEBUG" == true ]]; then
    echo "[DEBUG] Non-exited PIDs: ${non_exited_pids[*]}"
    echo
fi

# Find PIDs that exited but never had a start message (shouldn't happen in clean logs)
orphaned_exit_pids=()
if [[ "$DEBUG" == true ]]; then
    echo "[DEBUG] Finding PIDs that exited but never had a start message..."
fi
for pid in "${all_exited[@]}"; do
    if [[ "$DEBUG" == true ]]; then
        echo "  [DEBUG] Checking if PID $pid has a start log..."
    fi
    if [[ -z "${started_pids[$pid]+isset}" ]]; then
        if [[ "$DEBUG" == true ]]; then
            echo "    [DEBUG] PID $pid has NO start log, adding to orphaned_exit_pids"
        fi
        orphaned_exit_pids+=("$pid")
    else
        if [[ "$DEBUG" == true ]]; then
            echo "    [DEBUG] PID $pid HAS a start log"
        fi
    fi
done

if [[ "$DEBUG" == true ]]; then
    echo "[DEBUG] Orphaned exit PIDs: ${orphaned_exit_pids[*]}"
    echo
fi

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
    if [[ "$DEBUG" == true ]]; then
        echo "[DEBUG] Checking process status for PID $pid..."
    fi
    
    # Check if process exists and its state
    if ps -p "$pid" -o stat= 2>/dev/null 1>&2; then
        state=$(ps -p "$pid" -o stat= 2>/dev/null)
        # Remove leading/trailing whitespace
        state=$(echo "$state" | xargs)
        if [[ "$DEBUG" == true ]]; then
            echo "  [DEBUG] PID $pid exists with state: '$state'"
        fi
        
        if [[ "$state" == *"Z"* ]]; then
            if [[ "$DEBUG" == true ]]; then
                echo "    [DEBUG] PID $pid is a ZOMBIE"
            fi
            zombie_pids+=("$pid")
        else
            if [[ "$DEBUG" == true ]]; then
                echo "    [DEBUG] PID $pid is RUNNING"
            fi
            running_pids+=("$pid")
        fi
    else
        # Process doesn't exist
        if [[ "$DEBUG" == true ]]; then
            echo "  [DEBUG] PID $pid does NOT exist in process table (vanished)"
        fi
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
