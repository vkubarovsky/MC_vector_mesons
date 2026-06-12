# tmux Quick README

## Start / Attach

tmux new -s name        # start new session
tmux ls                 # list sessions
tmux attach -t name     # attach to session
tmux kill-session -t name

---

## Basic Controls (Prefix = Ctrl+b)

Ctrl+b d        # detach session
Ctrl+b c        # new window
Ctrl+b n        # next window
Ctrl+b p        # previous window
Ctrl+b w        # list windows

---

## Panes (Split Screen)

Ctrl+b %        # vertical split
Ctrl+b "        # horizontal split
Ctrl+b arrows   # move between panes
exit            # close pane

---

## Resize / Scroll

Ctrl+b Ctrl+arrow   # resize pane
Ctrl+b [            # scroll mode
q                   # exit scroll

---

## Rename

Ctrl+b ,        # rename window
Ctrl+b $        # rename session

---

## Minimal Workflow

tmux new -s work
# run your job
Ctrl+b d
# later:
tmux attach -t work

---

## Summary

- tmux keeps sessions alive after SSH disconnect
- allows multiple windows and splits
- essential for remote work
