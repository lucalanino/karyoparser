#!/usr/bin/env bash
# Fires after Edit/Write. If the edited file is in R/, reminds Claude to check docs.
f=$(jq -r '.tool_input.file_path // .tool_response.filePath // ""')
# Normalize backslashes to forward slashes (Windows paths)
f="${f//\\//}"

case "$f" in
  */R/*.R)
    jq -n --arg f "$f" '{
      "hookSpecificOutput": {
        "hookEventName": "PostToolUse",
        "additionalContext": "R source file edited: \($f)\nBefore finishing, verify documentation is up to date:\n- roxygen @param/@return for any changed exported function (check defaults, column lists)\n- CLAUDE.md architecture/API sections if internals or signatures changed\n- README.md if user-facing behavior changed\nRun devtools::document() if any roxygen changed."
      }
    }'
    ;;
esac
