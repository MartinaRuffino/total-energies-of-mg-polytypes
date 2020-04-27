(defun get-query (msg unit)
" Returns one of y, n or a, and queries if unit not a"
(interactive)
  (if (not (char-equal unit ?a))
     (progn (message msg) (setq unit (read-char)))
     (setq unit unit))
)

(defun f77-to-cray ()
"  Converts f77 programs to Cray"
   (interactive)
   (fortran-msm-intrsc-to-sngl (point-min) (point-max) nil)
)

(defun fortran-strip-long-lines ()
"   Strip columns beyond 72nd and trailing blanks in buffer."
   (interactive)
   (goto-char (point-min))
; Strip long lines ...
   (while (< (point) (point-max))
     (end-of-line) (if (> (line-len) 72) (delete-backward-char (- (line-len) 72)))
     (forward-char 1)
   )
; Strip trailing blanks ...
   (goto-char (point-min)) (replace-regexp "[ \t]+$" "")

  (kill-trailing-blank-lines)
)

(defun fortran-compare-call-args ()
"   Compare arguments of a subroutine call with subroutine args"
   (interactive)

; setup ...
   (beginning-of-line)
   (fortran-next-statement)
   (fortran-previous-statement)
   (setq thisbuf (current-buffer))
   (if (looking-at " +call +") (forward-sexp 1))(skip-chars-forward " \t")
   (setq subnam (buffer-substring (point) (progn (skip-chars-forward "a-z0-9") (point))))

; list all args of sub ...
   (save-excursion
     (find-tag subnam) (forward-sexp 2)
     (setq subrstrng (buffer-substring (+ 1 (point)) (progn (forward-sexp) (- (point) 1))))
     (pop-to-buffer "*compare*") (erase-buffer) (insert subrstrng "\n")
     (goto-char (point-min))
     (while (< (point) (point-max))
       (forward-sexp) (while (< (current-column) 20) (insert " "))
       (while (looking-at "[ \n\t,]") (delete-char 1)) (insert "\n")))

; list all args of call ...
     (setq subrstrng (buffer-substring (+ 1 (point)) (progn (forward-sexp) (- (point) 1))))
     (setq strngsiz 1)
     (pop-to-buffer "*compare*") (goto-char (point-min))
     (while (> strngsiz 0)
       (end-of-line) (setq thispt (point)) (insert subrstrng) (setq endpt (point))
       (goto-char thispt) (forward-sexp) (if (looking-at " *(") (forward-sexp))
       (while (and (looking-at "[ \t,]") (> endpt (point)))
	 (delete-char 1) (setq endpt (- endpt 1)))
       (if (and (looking-at "\n") (> endpt (point)))
	   (progn (delete-char 7) (setq endpt (- endpt 7))
		  (while (and (looking-at "[ \t,]") (> endpt (point)))
		    (delete-char 1) (setq endpt (- endpt 1)))))
       (setq subrstrng (buffer-substring (point) endpt))
       (setq strngsiz (- endpt (point)))
       (delete-char (- endpt (point)))
; (insert-string strngsiz)  (insert-string subrstrng)
       (next-line 1)

; (setq strngsiz 0)
     )
)

(defun kill-trailing-blank-lines ()
"  Kills blanks, tabs and \n at end of file, with last char \n"

   (interactive)

   (goto-char (point-max))
   (skip-chars-backward " \t\n")
   (delete-region (point) (point-max))
   (insert "\n")
)

(defun fortran-file-as-subnam ()
"  sets the file name to first subnam, extension .f"

   (interactive)

   (re-search-forward "^ +\\(integer\\|double precision\\|logical\\|complex\\|complex\\*16\\|real\\|real\\*8\\)* +function\\|^ +subroutine\\|^ +program")
   (skip-chars-forward " \t")
   (setq nwnam (downcase (buffer-substring (point)
				 (progn (skip-chars-forward "a-zA-Z0-9")
					(point)))))
   (set-visited-file-name (concat nwnam ".f"))
)

(defun line-len ()
"  Returns length of line point is on."
   (interactive)
   (- (save-excursion (end-of-line) (point)) (save-excursion (beginning-of-line) (point))))

(defun cms-to-cray ()
"  Converts MSM's CMS LMTO programs to MSM's Cray LMTO programs."
   (interactive)
   (mark-whole-buffer)
   (replace-regexp "\\b\\$" "J7")
   (mark-whole-buffer)
   (while (< (point) (point-max))
     (while (looking-at "^[Cc*]") (next-line 1))
; cast conversions
     (kill-line 1) (yank) (previous-line 1)
     (if (/= 0
         (fortran-msm-intrsc-to-sngl (point) (save-excursion (end-of-line) (point)) nil))
	 (progn (beginning-of-line) (yank) (previous-line 1)
		(insert-string "C>") (delete-char 2)))
     (if (looking-at ".+INTEGER (J7)") (progn
        (kill-line 1) (yank) (yank)
        (previous-line 2) (insert-string "C>") (delete-char 2)
        (next-line 1) (re-search-forward "(J7)") (replace-match "(I)")
     ))
     (if (looking-at " +CALL ERRSET") (progn
	(insert-string "C>") (delete-char 2)
     ))
     (beginning-of-line) (next-line 1)
  )
)

(defun fortran-msm-intrsc-to-sngl (startposn posn que)
  "converts intrinsic fortran functions to sngl in region (STARTPOSN, POSN).
   If QUE (third arg) is t, asks interactively.
   Returns change in size of region."
  (interactive)
; Set up for query replace
   (if que (setq query ??) (setq query ?a))
   (setq savpos posn)

; Handle all cases where dfnct -> fnct
   (goto-char startposn)
   (while (re-search-forward "\\bd\\(cmplx\\|conjg\\|real\\)\\b"  posn t) (progn
      (setq query (get-query "Convert?" query))
      (if (or (char-equal query ?a) (char-equal query ?y))
         (progn (backward-word 1) (delete-char 1)
            (setq posn (- posn 1))
   ))))

; Handle all cases where cdfnct -> cfnct
    (goto-char startposn)
    (while (re-search-forward
            "\\bcd\\(exp\\|abs\\)\\b"  posn t) (progn
      (setq query (get-query "Convert?" query))
      (if (or (char-equal query ?a) (char-equal query ?y))
         (progn (backward-word 1) (forward-char 1) (delete-char 1)
              (setq posn (- posn 1))
   ))))

; Handle all cases where dfnct -> afnct
   (goto-char startposn)
   (while (re-search-forward
           "\\bd\\(imag\\)\\b" posn t) (progn
      (setq query (get-query "Convert?" query))
      (if (or (char-equal query ?a) (char-equal query ?y))
         (progn (backward-word 1) (delete-char 1) (insert-string "a"))
   )))

 (- posn savpos))                    ; Return with change in region size

(defun cray-to-cms ()
"  Converts MSM's CRAY LMTO programs to MSM's CMS LMTO programs."
   (interactive)
   (goto-char (point-min))
   (replace-regexp "\\bJ7" "$")
   (goto-char (point-min))
   (while (< (point) (point-max))
     (if (looking-at "C> +CALL ERRSET") (progn
	(insert-string "  ") (delete-char 2)
     ))
     (if (looking-at "^[Cc*][>]") (progn 
	(insert-string "  ") (delete-char 2)
        (next-line 1) (beginning-of-line)
	(fortran-intrsc-to-double (point) (save-excursion (end-of-line) (point)) nil)
        (beginning-of-line)
	(upcase-region (point) (save-excursion (end-of-line) (point)))
        (previous-line 1)
	(if (not (looking-at "\\(.+\\)$\n\\1"))
	    (message "----- unmatched conversion line ----"))
        (next-line 1) (kill-line 1) (previous-line 1)
     ))
     (beginning-of-line) (next-line 1)
  )
)

(defun cms-main-to-apollo ()
" fix apollo bug for too large local stack"
   (interactive)

; Check to convert to unix fmain
     (goto-char (point-min))
     (insert-string "      subroutine fmain\n")

     (if (re-search-forward "^ +PROGRAM" (point-max) t)
	 (progn (beginning-of-line)
		(insert-string "C>") (delete-char 2)))
     (goto-char (point-min))

; Put ERRSET as first procedure call
     (if (re-search-forward "^ +CALL ERRSET" (point-max) t) (progn
; old: comment it out...	(insert-string "C>") (delete-char 2)
	(beginning-of-line)
	(kill-line 1)
	(fortran-first-executable-statement)
	(yank)
	(setq erset nil)
     ))
     (goto-char (point-min))

  (if (re-search-forward "\\(^ +common */w/\\)" (point-max) t)
      (progn
	(beginning-of-line) (next-line 1)
	(insert-string "      common /static/\n")
	))

  (if (re-search-forward "\\(^ +CALL WKINFO\\)" (point-max) t)
      (progn
	(beginning-of-line)
	(delete-char 1) (insert-string "C")
	(beginning-of-line) (next-line 1)
	(insert-string "      call fexit(0,' ')\n")))
)

(defun cms-to-apollo ()
"  Converts MSM's CMS LMTO programs to Apollo LMTO programs."
   (interactive)
   (goto-char (point-min))
   (fortran-first-executable-statement)

   (setq erset t)
   (while (< (point) (point-max))
     (while (looking-at "^[Cc*]") (fortran-next-statement))

; Check to convert to unix fmain
;     (if (looking-at "^ .+PROGRAM")
;	 (save-excursion 
;	   (delete-char 1) (insert-string "C")
;	   (goto-char (point-min))
;	   (insert-string "      subroutine fmain\n")
;	   (fortran-first-executable-statement)
;	   (beginning-of-line)
;	   (insert-string "      save\n")))


; File opening for control file
     (if (looking-at ".+CALL IOCF") (save-excursion (fortran-add-fopn "10" t t)))

; File opening for chd files
     (if (looking-at ".+CALL IOCHD") (progn
        (re-search-forward "chd.?(")
	(backward-char 1) (forward-sexp 1) (backward-word 1)
	(setq unit (buffer-substring (point)
				     (save-excursion (re-search-forward "[0-9]+") (point))))
	(save-excursion (fortran-add-fopn unit nil t))
	))

; File opening for atomic file
     (if (looking-at ".+CALL AIOGEN") (progn
        (re-search-forward "gen(")
	(setq unit (buffer-substring (point) (save-excursion (re-search-forward "[a-z0-9()]+") (point))))
	(beginning-of-line)
	(save-excursion (fortran-add-fopn unit nil nil))
        (insert-string "      ifi = fopna(clabel(ic),ic+10,0)\n")
	(fortran-next-statement)
        (save-excursion 
        (if (re-search-backward "^ +do +[0-9]+ +ic" (point-min) t) (fortran-next-block))
        (insert-string "      call fclose(ifi)\n"))
     ))

; File opening for dpdump ...
     (if (looking-at " +call dpdump") (save-excursion 
       (end-of-line) (backward-word 1)
	(setq unit (buffer-substring (point) (save-excursion (re-search-forward "[0-9]+") (point))))
	(fortran-add-fopn unit nil t)
     ))

; File opening for rewind ...
;     (if (looking-at " +rewind") (save-excursion 
;       (end-of-line) (backward-word 1)
;	(setq unit (buffer-substring (point) (save-excursion (re-search-forward "[0-9]+") (point))))
;	(fortran-add-fopn unit nil t)
;     ))

; File opening for other files
     (if (looking-at " .+\\(read *( *[^5*]\\\|write *( *[^-6*]\\)") (save-excursion 
	(re-search-forward " .+\\(read\\\|write\\) *( *")        
	(setq unit (buffer-substring (point) (save-excursion (re-search-forward "[0-9]+") (point))))
	(fortran-add-fopn unit (looking-at "[0-9]+ *, *[1-9]") t)
     ))

     (beginning-of-line) (fortran-next-statement)
  )
  (setq erset t)
)

(defun fortran-add-fopn (unit ifformat ifopn)
" Adds fopn call to current function, opening UNIT, formatted.
  Second arg is true if file is formatted, false if unformatted
  Third argument is nil if only declare integer function fopna.' "
   (interactive)
   (fortran-first-executable-statement)

; First, declare fopna as integer
   (fortran-previous-statement)
   (while (looking-at " +data") (fortran-previous-statement))
   (if (looking-at " +integer fopna") nil
; Otherwise, look for commented-out declaration ...
      (if (looking-at "C# +integer fopna")
	  (progn (insert-string "  ") (delete-char 2))
; Otherwise, insert one ...
	  (progn (fortran-next-statement) (insert-string "      integer fopna\n")))
   )

; Next, add "call fopna"
   (if (and ifopn (cms-file-name (string-to-int unit))) (progn
     (fortran-first-executable-statement)
    (while (and (looking-at " +o = fopna(")
                (save-excursion (re-search-forward "fopna([^,]+,")
		   (string< (buffer-substring (point) (progn (skip-chars-forward "0-9") (point))) unit)))
       (next-line 1)
     )
;    Insert the call
     (insert-string (concat "      o = fopna('" (cms-file-name (string-to-int unit)) "',"  unit "," (if ifformat 0 4) ")\n"))
     (previous-line 1)
   
; Check and eliminate duplication ...
     (if (looking-at "\\(.+\\)$\n\\1") (kill-line 1))
  ))
)


(defun cms-file-name (arg)
"  Converts file number into name, MSM's conventions."
   (interactive)

  (cond
        ((= arg  9) "rfit")
        ((= arg 10) "lc")
        ((= arg 60) "hdt")
        ((= arg 61) "hfit")
        ((= arg 62) "hus")
        ((= arg 63) "hsxi")
        ((= arg 64) "hint")
        ((= arg 65) "hnt0")
        ((= arg 66) "hst0")
        ((= arg 67) "hst1")
        ((= arg 68) "hst2")
        ((= arg 69) "hpvs")
        ((= arg 71) "log")
        ((= arg 72) "rho")
        ((= arg 73) "pert")
        ((= arg 74) "rhoi")
        ((= arg 75) "poti")
        ((= arg 76) "dos")
        ((= arg 77) "qpts")
        ((= arg 78) "rhoe")
        ((= arg 79) "rh")
        ((= arg 80) "rl")
        ((= arg 81) "pldt")
        ((= arg 88) "ms")
        ((= arg 89) "sv")
        ((= arg 91) "wrk1")
        ((= arg 92) "wrk2")
        ((= arg 93) "syml")
        ((= arg 94) "bnds")
        ((= arg 95) "rfmt")
  )
)
