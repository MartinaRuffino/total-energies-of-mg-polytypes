; --- macro fortran-call-tree ---
(defun fortran-call-tree (&optional nest outmode)
   " Tree for subroutine call.
     Optional NEST sets upper bound to the nesting of calls.
     From a program, outmode nonnil lists files procedures are in.
     A * after name means tag not in table; ... means name already found."

    (interactive "P")
    (if (not nest) (setq nest 7))

; --- Setup ---
   (progn
   (fortran-mode) (visit-tags-table "./")
   (setq buffer (car (buffer-list))) (setq bnam (buffer-name buffer))
   (pop-to-buffer "*Tree*") (kill-buffer "*Tree*")
   (pop-to-buffer "*Tree*") (pop-to-buffer bnam)
   )

   (fortran-call-tree2 bnam nest outmode "  "
     (point) (save-excursion (end-of-fortran-subprogram) (point)))

   (pop-to-buffer "*Tree*") (pop-to-buffer bnam)

 )

(defun fortran-call-tree2 (bnam level outmode spaces pointnow eos)
" Called by fortran-call-tree.
  bnam       : buffer name
  level      : number of nesting levels remaining
  outmode    : How messaging is handled?
  spaces     : spacing to precede text in tree output
  pointnow   : current point
  eos        : position marking routine end
  skip-prior : If tag already parsed,  ?
"

  (interactive "P")

  (let ((tag-already-found nil) tagnam strn ncallarg nsubarg cln)

; --- For each call do ---
  (while (re-search-forward "^[ 0-9]....+\\bcall +\\([A-Za-z0-9_]+\\) *(" eos t)
    ;; Preparatory step: see whether call has tag, and whether already found
    (progn
      ; tag
      (setq tagnam (buffer-substring-no-properties (match-beginning 1) (match-end 1)))
      (backward-char 1)
      (setq ncallarg (fortran-count-call-args eos))
      (setq pointnow (point))

      ;; File containing tag
      (setq tagfil (if (is-tag tagnam)
		       (progn (visit-tags-table-buffer) (file-of-tag)) "*"))

      ;; Determine whether tag has already been found; make initial string for *Tree*
;      (setq tagfil "*")
      (if (string= tagfil "*")
	  (setq tag-already-found nil) ; no tag; nothing to do

	;; search in *Tree* for prior occurence of tag
	(set-buffer "*Tree*") (goto-char (point-max))
	(setq tag-already-found
	      (save-excursion (re-search-backward
			       (concat "\\b" tagnam "\\b") (point-min) t)))
	)

      ;; First cut at printout for *Tree*
      (setq strn (concat tagnam "  " (if tag-already-found " ... " tagfil)))

      ;; return to buffer
      (set-buffer bnam) (setq pointnow (point))

      ) ; end of preparatory step

;    (snit)

    ;; Output info for this call; possibly look for nesting of tree
    (if (and (is-tag tagnam) (> level 0))
	;; Branch where tag has found.
	;; Follow tag to check arg count; recursively follow tree
	(progn
	  ;; current line number for future reference
	  (setq cln (count-lines (point-min) (point)))
	  (if (= (current-column) 0) (setq cln (+ cln 1)))

	  ;; Go to file containing tag
	  (find-tag tagnam)

	  ;; buffer name of tag file
	  (setq buffer (car (buffer-list)))  (setq tnam (buffer-name buffer))

	  ;; Count arguments associated with tag
	  (if (not (re-search-forward "^ .+\\b\\(subroutine\\|function\\|entry\\) +\\([A-Za-z0-9_]+\\) *("
			     (save-excursion (fortran-next-statement) (point)) t))
	      (setq nsubarg 0)
	    (backward-char 1)
	    (setq nsubarg (fortran-count-call-args (save-excursion (forward-sexp 1) (point)))))

	  ;; print number of args, or if mismatch, also mismatch info
	  (if (= ncallarg nsubarg)
	      (setq strn (concat strn "  " (int-to-string nsubarg))) ; append number of args
					; append number of args, line number of mismatch
;	    (snit)
	    (setq strn (concat strn "  " (int-to-string nsubarg) " mismatch " bnam " " (int-to-string cln)))
	    )

	  ;; Append information to *Tree* buffer
	  (set-buffer "*Tree*") (goto-char (point-max))
	  (if (not outmode) (message "%s%s" spaces strn))
;	  (if outmode (message "%s" (buffer-file-name)))
	  (insert spaces strn "\n")

;	  (snit)

	  ;; recursively continue at next nesting level
	  (if tag-already-found nil
	    (set-buffer buffer)
	    (fortran-call-tree2 tnam (- level 1) outmode (concat "  " spaces) (point)
				(save-excursion (end-of-fortran-subprogram) (point)))
	    )
	  ) ; end of branch where tag is available

      ;; Branch where no tag to follow: add to *Tree* without checking arg count
      (set-buffer "*Tree*") (goto-char (point-max))
      (insert spaces strn "\n")
      (if (not outmode) (message "%s%s" spaces strn))
      (if outmode (message "%s" (buffer-file-name)))
      ) ; end of all branches handling output for this subroutine call

      ;; Continue in current buffer
      (pop-to-buffer bnam) (goto-char pointnow)

    ) ; search for subroutine calls exhausted
))

(defun fortran-count-call-args (maxpt)
" Counts number of arguments enclosed by next (..)
  Language assumed to be fortran."

    (interactive "P")

;    (snit)
    (let (pointnow first last (nargs 0))
      (setq pointnow (point))
;      (snit)
      (if (not (search-forward "(" maxpt)) nil
	(setq first (point)) (backward-char 1)
	(setq last (progn (forward-sexp 1) (backward-char 1) (point)))
	(goto-char first)
	(while (< (point) last)
	  (setq nargs (+ nargs 1))
	  (fortran-forward-expr 1) (forward-char 1) (skip-chars-forward " "))
	)
    (goto-char pointnow)
    nargs
))
