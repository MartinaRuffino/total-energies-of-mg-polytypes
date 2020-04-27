; --- macro fortran-list-subroutine-args (21 Nov 13) ---
(defun fortran-list-subroutine-args (&optional arg)
"   Document arguments of a subroutine.
    Scalars in file vars are used where available.
    Structures elements are identified grouped according to whether they are
    read, stored, allocated (in the case of pointers), or passed as arguments to subprograms.

    Optional ARG non-nil causes a (recursive) search to find structure elements in routines
    called by this one.
"
   (interactive "P")

   (let (thisbuf strngsiz thispt pparpt intpt dppnt endpt subrstrng strxnam subnam substart subend founds strxelt ifstore
	 sizcol r start end thisarg varfilnam varbuf switchfilnam (strxline 0)
	 (prefix-is-C-u
 		  (and arg (not (equal arg (prefix-numeric-value arg))))))

; .. Hang onto buffer name
     (setq thisbuf (current-buffer))
; .. file containing default definitions for variables
     (setq varfilnam "~/lm/misc/vars")
; .. file containing default definitions for switches
     (setq switchfilnam "~/lm/misc/command-line-switches")
; .. Find start of subprogram
     (beginning-of-line)
     (fortran-next-statement)
     (fortran-previous-statement)
     (if (not (looking-at " +\\(subroutine\\|integer function\\)"))
	 (beginning-of-fortran-subprogram))
     (setq substart (point))
; .. Find end of subprogram
     (save-excursion
       (end-of-fortran-subprogram)
       (setq subend (point)))

; .. Pick up argument list
     (save-excursion (setq subrstrng (re-search-forward "(" subend 0)))
     (save-excursion
       (setq subrstrng (buffer-substring-no-properties
			(progn (re-search-forward "(" subend 0)
			       (backward-char 1) (+ (point) 1))
			(progn (forward-sexp) (- (point) 1)))))
     (setq substrng subrstrng)

; .. Build up command line arguments in *cmdlst*
     (pop-to-buffer "*cmdlst*") (erase-buffer)

; .. Build up comments in *arglst*
     (pop-to-buffer "*arglst*") (erase-buffer) (insert subrstrng "\n")
					; Clean up around newlines
     (setq sizcol 0)
     (goto-char (point-min))
     (while (progn (forward-line 1) (< (point) (point-max)))
	(delete-char 6)
	(while (looking-at "[ \t]") (delete-char 1))
	(insert ". "))
					; Put arguments into columns; get max column size
     (goto-char (point-min))
     (while (< (point) (point-max))
       (forward-sexp) (setq sizcol (max sizcol (current-column)))
       (while (looking-at "[ \n\t,]") (delete-char 1)) (insert "\n"))
     (setq sizcol (+ sizcol 1))

					; Replace "^ \." with a newline now
     (goto-char (point-min))
     (while (< (point) (- (point-max) 0))
       (if (looking-at "^\\.") (progn (delete-char 2) (insert "\n")))
       (next-line 1))

; ...For each structure, collect inputs and outputs
     (goto-char (point-min))
     (while (< (point) (- (point-max) 0))

;       (snit)
       ; Move to next line starting with s_...
       (if (not (looking-at "s_"))
	   (forward-line 1)
	 (setq subrstrng (buffer-substring-no-properties
			  (progn (beginning-of-line) (+ (point) 2))
			  (progn (end-of-line) (point))))
	 (setq strxnam (concat "s_" subrstrng))
	 (pop-to-buffer thisbuf)  (goto-char substart)

;        (snit)
         ; If this struc is passed in a sub call, set up comments
	 (setq founds nil)
	 (if (not (or
	       (re-search-forward (concat strxnam ".*%") subend t)
	       (re-search-forward (concat strxnam " *[,)]" ) subend t)
	       )) nil
	   ; case some element of structure used by this routine
	   (pop-to-buffer "*arglst*")
	   (progn (setq founds (point))
		  (end-of-line)
		  (insert "\n  Elts read: \n  Stored:    \n  Allocated: \n  Elts passed:\n  Passed to: "))

	   (pop-to-buffer thisbuf)  (goto-char substart))
;	   (pop-to-buffer thisbuf) (goto-char (point-min)))

         ;; ... For each each occurence where an element is packed or unpacked, add to declarations
;	 (snit)
	 (fortran-find-strux-elements thisbuf founds varfilnam strxnam substart subend arg)

	 ;; Clean-up structure descriptions
;	 (snit)
	 (pop-to-buffer "*arglst*")
	 (if (not founds)
	     (progn
	       (insert " (not used)")
	       (forward-line 1))  ; move to next argument

					; Eliminate duplicates in 'Stored:' line
	   (goto-char founds)
	   (search-forward "Stored:")
	   (if (looking-at " *$") nil
	     (while
		 (re-search-forward "\\( [^ ]+\\)\\b.*\\(\\1\\)\\b" (save-excursion (end-of-line) (point)) t)
	       (replace-match "" nil t nil 2) (beginning-of-line)))

					; Eliminate duplicates in 'Allocated:' line
	   (goto-char founds)
	   (search-forward "Allocated:")
	   (if (looking-at " *$") nil
	     (while
		 (re-search-forward "\\( [^ ]+\\)\\b.*\\(\\1\\)\\b" (save-excursion (end-of-line) (point)) t)
	       (replace-match "" nil t nil 2) (beginning-of-line)))

					; Eliminate duplicates in 'Elts passed:' line
	   (goto-char founds)
	   (search-forward "Elts passed:")
	   (if (looking-at " *$") nil
	     (while
		 (re-search-forward "\\( [^ ]+\\)\\b.*\\(\\1\\)\\b" (save-excursion (end-of-line) (point)) t)
	       (replace-match "" nil t nil 2) (beginning-of-line)))

					; Eliminate duplicates in 'Passed to:' line
	   (goto-char founds)
	   (search-forward "Passed to:")
	   (if (looking-at " *$") nil
	     (while
		 (re-search-forward "\\( [^ ]+\\)\\b.*\\(\\1\\)\\b" (save-excursion (end-of-line) (point)) t)
	       (replace-match "" nil t nil 2) (beginning-of-line)))

				; Eliminate duplicates in 'Read:' line
	   (goto-char founds)
	   (search-forward "Read:")
	   (if (looking-at " *$") nil
	     (setq strxelt)
	     (while
		 (re-search-forward "\\( [^ ]+\\)\\b.*\\(\\1\\)\\b" (save-excursion (end-of-line) (point)) t)
	       (setq strxelt (concat strxelt
				     (buffer-substring-no-properties (match-beginning 2) (match-end 2))))
	       (replace-match "" nil t nil 2)
	       (beginning-of-line))
					; Write out duplicates for information
;; 	    (if (not strxelt) nil
;; 	      (search-forward "Stored:")
;; 	      (beginning-of-line) (insert "  Duplicate:" strxelt "\n"))

	     )

;	 (snit)

					; Shorten 'Read:' line to 67 characters
	  (goto-char founds)
	  (search-forward "Read:")
	  (while (> (save-excursion (end-of-line) (current-column)) 67)
	    (beginning-of-line) (forward-char 66) (skip-chars-forward " ") (if (looking-at "\\B") (backward-word 1))
	    (insert "\n              ")
	    )
	  (beginning-of-line) (if (looking-at "  Elts read: +$")
				  (save-excursion (end-of-line) (insert " *")))

					; Shorten 'Elts passed:' line to 67 characters
	  (goto-char founds)
	  (search-forward "Elts passed:")
	  (while (> (save-excursion (end-of-line) (current-column)) 67)
	    (beginning-of-line) (forward-char 66) (skip-chars-forward " ") (if (looking-at "\\B") (backward-word 1))
	    (insert "\n              ")
	    )
	  (beginning-of-line) (if (looking-at "  Elts passed:$")
				  (save-excursion (end-of-line) (insert "*"))
				(if (looking-at "  Elts passed: ") (progn (forward-char 14) (delete-char 1))))

					; Shorten 'Allocated to:' line to 67 characters
	  (goto-char founds)
	  (search-forward "Allocated:")
	  (while (> (save-excursion (end-of-line) (current-column)) 67)
	    (beginning-of-line) (forward-char 66) (skip-chars-forward " ") (if (looking-at "\\B") (backward-word 1))
	    (insert "\n              ")
	    )
	  (beginning-of-line) (if (looking-at "  Allocated: +$")
				  (save-excursion (end-of-line) (insert " *")))

					; Shorten 'Stored:' line to 67 characters
	  (goto-char founds)
	  (search-forward "Stored:")
	  (while (> (save-excursion (end-of-line) (current-column)) 67)
	    (beginning-of-line) (forward-char 66) (skip-chars-forward " ") (if (looking-at "\\B") (backward-word 1))
	    (insert "\n              ")
	    )
; 	  (snit)
	  (beginning-of-line) (if (looking-at "  Stored: +$")
				  (save-excursion (end-of-line) (insert " *")))

					; Shorten 'Passed to:' line to 67 characters
	  (goto-char founds)
	  (search-forward "Passed to:")
	  (while (> (save-excursion (end-of-line) (current-column)) 67)
	    (beginning-of-line) (forward-char 66) (skip-chars-forward " ") (if (looking-at "\\B") (backward-word 1))
	    (insert "\n              "))

	  (beginning-of-line) (if (looking-at "  Passed to: +$")
				  (save-excursion (end-of-line) (insert " *")))
	  (beginning-of-line) (if (looking-at "  Passed to:") (next-line 1))
	  ) ; end of cleanup.  Current buffer should be "*arglst*"
	 (beginning-of-line)))
     (delete-blank-lines)

; .. Open file 'vars'
     (find-file varfilnam)
     (setq varbuf (current-buffer))

; .. Find description for any names that match
     (pop-to-buffer "*arglst*")
     (goto-char (point-min))
     (while (< (point) (- (point-max) 0))
       (if (or (looking-at "^$") (looking-at "^ ")) (next-line 1)
       (setq subrstrng (buffer-substring-no-properties
			(progn (beginning-of-line) (point))
			(progn (end-of-line) (point))))
       (pop-to-buffer varbuf)
       (goto-char (point-min))
; (snit)
       (if (re-search-forward (concat "^\\(w(\\)?" subrstrng "\\b[^_]") nil t)
	   (progn (beginning-of-line)
		  (while (re-search-forward (concat "^\\(w(\\)?" subrstrng "\\b[^_]") nil t)
		  (beginning-of-line)
		  (setq strxnam (buffer-substring-no-properties
				   (progn (beginning-of-line) (point))
				   (progn (end-of-line) (point))))
		  (pop-to-buffer "*arglst*")
		  (beginning-of-line) (setq thispt (point))
		  (insert strxnam "\n")
					;this segment formats inserted text
		  (while
		      (progn
			(setq endpt (point))
			(goto-char thispt)
			(re-search-forward "%N" endpt t))
		    (delete-char -2) (insert "\n       ") (goto-char (+ 6 endpt))
		    t) (goto-char endpt)
		  (pop-to-buffer varbuf))
		  (pop-to-buffer "*arglst*")
		  (kill-line 1))
	 (progn (pop-to-buffer "*arglst*") (next-line 1) (beginning-of-line)))))

; ... Convert declarations into comment lines
     (goto-char (point-min))
     (insert "C-\n")
     (insert "C ----------------------------------------------------------------------\n")
     (insert "Ci Inputs\n")
     (while (< (point) (- (point-max) 0))
       (if (looking-at "^$") (progn (insert "C") (next-line 1) (beginning-of-line)))
       (if (looking-at "^s_")
	   (progn
	     (if (= strxline 0)  (insert "Cio Structures\n"))
	     (setq strxline 4)
;	     (snit)
	     (insert "Cio  ") (next-line 1) (beginning-of-line))
	 (if (and (> strxline 1) (not (looking-at "^ "))) (setq strxline 1))
	 (if (and (> strxline 1) (looking-at "^  Elts read")) (setq strxline 3))
	 (if (and (> strxline 1) (looking-at "^  Stored")) (setq strxline 2))
	 (if (and (> strxline 1) (looking-at "^  Allocated")) (setq strxline 2))
	 (if (and (> strxline 1) (looking-at "^  Argument")) (setq strxline 5))
	 (if (and (> strxline 1) (looking-at "^  Elts passed")) (setq strxline 5))
	 (progn
	   (if (= strxline 5) (insert "Cio  ")
	     (if (= strxline 3) (insert "Ci   ")
	       (if (= strxline 2) (insert "Co   ")
		 (insert "Ci   ")))))
		(next-line 1) (beginning-of-line))
	 )
     (insert "Co Outputs\nCs Command-line switches\n")

; ... Pick up command-line arguments, insert as comment lines
;    (snit)
     (pop-to-buffer thisbuf)  (goto-char substart)
     (fortran-find-cmdargs thisbuf varfilnam substart subend arg)
     (pop-to-buffer "*cmdlst*")
     (goto-char (point-min))
;     (snit)
     (if (looking-at " *$") (pop-to-buffer "*arglst*")
       (insert " ") (backward-char 1)
       (while (re-search-forward " \\([^ ]+ \\).*\\(\\1\\)" (save-excursion (end-of-line) (point)) t)
	 (progn
	   (replace-match "" nil t nil 2) (beginning-of-line))
	 )
					; ... Put arguments into columns
       (beginning-of-line) (delete-char 1)
       (replace-string " " "\n")
       (sort-lines nil (point-min) (point-max))

					;  Setup for finding descriptions of switches
       (find-file switchfilnam)
       (setq varbuf (current-buffer))
       (pop-to-buffer "*cmdlst*")
					; ... Max column size
       (goto-char (point-min))
       (setq sizcol 0)
       (while (< (point) (point-max))
	 (end-of-line) (setq sizcol (max sizcol (current-column))) (forward-char 1))
       (setq sizcol (+ sizcol 1))

					; Embellish --switch
;       (snit)
       (goto-char (point-min))
       (while (< (point) (point-max))
	 (setq subrstrng (buffer-substring-no-properties
			  (progn (beginning-of-line) (point))
			  (progn (end-of-line) (point))))
	 (end-of-line)
	 (while (< (current-column) sizcol) (insert " "))
	 (insert ":") 		; ... Put ":" at end of line
	 (pop-to-buffer varbuf)  ;Check for default definition for this switch
	 (goto-char (point-min))
	 (if (not (re-search-forward subrstrng nil t)) (pop-to-buffer "*cmdlst*")
	   (skip-chars-forward "^[:]")
	   (forward-char 1)
	   (setq subrstrng (buffer-substring-no-properties
			    (point) (progn (end-of-line) (point))))
	   (pop-to-buffer "*cmdlst*")
	   (insert subrstrng) 		; Append description
	   )
	 (forward-char 1)
	 )
					; translate into fortran comments
       (goto-char (point-min))
       (while (< (point) (- (point-max) 0))
	 (progn (insert "Cs   ") (next-line 1) (beginning-of-line)))
       (setq subnam (buffer-substring-no-properties (point-min) (point-max)))
       (pop-to-buffer "*arglst*")
       (insert subnam)
       )

;     (snit)

; ... Remainder of header
     (insert "Cl Local variables\nCl         :\nCr Remarks\nCr   \nCu Updates\n")
     (insert "Cu   " (add-doc-iso8601-time-string) " \n")
     (insert "C ----------------------------------------------------------------------\n")
     (setq pparpt (point))
     (insert "      implicit none\nC ... Passed parameters\n" "      integer," substrng "\n")

					; declare as double [a-hp-z]
     (goto-char pparpt)
     (search-forward "integer")
     (if (looking-at " *$") nil
					;for now, condense into a single line
       (save-excursion (while (re-search-forward "\n     . *" nil t) (replace-match "" nil t)))

       (setq strxelt)
       (while
	   (save-excursion (re-search-forward "\\(,[a-hA-Hp-zP-Z][^,]*\\)\\b" (save-excursion (end-of-line) (point)) t))
	 (setq strxelt (concat strxelt
			       (buffer-substring-no-properties (match-beginning 1) (match-end 1))))
	 (replace-match "" nil t nil 1))

					; replace "integer," with "integer "
       (goto-char pparpt) (search-forward "integer") (delete-char 1) (insert " ")

					; Write out double precision variables
       (if (not strxelt) nil
	 (goto-char pparpt) (search-forward "integer") (forward-line 1)
	 (setq endpt (point))
	 (insert "      double precision" strxelt "\n"))

					; Shorten integer line to 72 characters
       (goto-char pparpt) (search-forward "integer")
       (while (> (save-excursion (end-of-line) (current-column)) 72)
	 (beginning-of-line) (forward-char 71) (skip-chars-forward ",") (if (looking-at "\\B") (backward-word 1))
	 (insert "\n     .        "))

					; Pull out arguments that are structures
       (goto-char pparpt) (setq strxelt)
       (if (not (search-forward "double precision" nil t)) nil
	 (while (re-search-forward "\\(,[a-hA-Hp-zP-Z][^,]*\\)\\b" (save-excursion (end-of-line) (point)) t)
	   (setq strxnam (buffer-substring-no-properties (match-beginning 1) (match-end 1)))
	   (if (not (save-excursion (goto-char (point-min)) (re-search-forward (concat (substring strxnam 1) ".+\nCi     Elts read") nil t))) nil
	     (search-backward strxnam) (replace-match "" nil t)
	     (setq strxelt      (concat strxelt "      type(str_" strxnam "):: " strxnam "\n"))))
	 (if (not strxelt) nil
	   (forward-line 1)
	   (insert strxelt "\n")))

					; Shorten double-precison line to 72 characters
       (goto-char pparpt)
       (while (search-forward "double precision" nil t)
	 (if (looking-at "$") nil
	   (delete-char 1) (insert " ")
	   (while (> (save-excursion (end-of-line) (current-column)) 72)
	     (beginning-of-line) (forward-char 71) (skip-chars-forward ",") (if (looking-at "\\B") (skip-chars-backward "a-zA-Z0-9_()"))
	     (insert "\n     .                 "))))

					; Clean up structure declarations
       (goto-char pparpt)
       (if (re-search-forward "\\(type(str_\\)\\(,s_\\)\\(.+\\),s" nil t)
	   (progn (beginning-of-line)
		  (insert "\nC ... For structures\n")
		  (insert "      include 'structures.h'\n")))

       (goto-char pparpt)
       (while (re-search-forward "\\(type(str_\\)\\(,s_\\)\\(.+\\),s" nil t)
	 (replace-match "\\1\\3s"))

       (goto-char pparpt)
       (if (re-search-forward "type(str_spec):: s_spec" nil t)
	   (insert "()"))

       (goto-char pparpt)
       (if (re-search-forward "type(str_site):: s_site" nil t)
	   (insert "()"))

       (goto-char (point-max)))

     (setq intpnt (point))


; --- List subroutines called by this routine (external declarations) ---
     (pop-to-buffer thisbuf)
     (setq strxelt)
     (goto-char (point-min))
     (while (re-search-forward "^[ \t].+[ \t]call +\\([a-zA-Z0-9_]+\\)" nil t)
       (setq strxelt (concat strxelt " " (buffer-substring-no-properties (match-beginning 1) (match-end 1)))))
     (pop-to-buffer "*arglst*")

;    Poke external calls into *arglst*
     (if (not strxelt) nil
       (insert "C ... External calls\n      external")
       (setq thispt (point))
       (save-excursion (insert strxelt "\n"))
;      replace " " with ","
       (goto-char thispt)
       (while (search-forward " " nil t) (replace-match "," nil t))

;      Eliminate duplicate names in line
       (goto-char thispt)
					; Sort entries in line
       (sort-regexp-fields nil "[^,]+" "\\&" thispt (save-excursion (end-of-line) (point)))
       (goto-char thispt)
       (while
	   (re-search-forward ",\\([^,]+\\)\\(,\\1\\)" (save-excursion (end-of-line) (point)) t)
	 (replace-match "" nil t nil 2) (goto-char thispt))
       (goto-char thispt) (delete-char 1) (insert " ")
					; Shorten line to 72 characters
       (while (> (save-excursion (end-of-line) (current-column)) 72)
	 (beginning-of-line) (forward-char 71) (skip-chars-forward ",") (if (looking-at "\\B") (backward-word 1))
	 (insert "\n     .         "))
       (forward-line 1))

; --- Local Variable declarations ---
     (goto-char intpnt)
     (insert "C ... Local parameters\n      integer ")
     (setq intpnt (point))
     (insert "\n      double precision ")
     (setq dppnt (point))
     (insert "\n")
;      (pop-to-buffer thisbuf) (goto-char (point-min))
;      (while (re-search-forward "\\b\\([a-zA-Z0-9_]+\\)\\b" nil t)
; 	   (snot)
;        (if (not (save-excursion (beginning-of-line) (looking-at " \t"))) nil
; 	 ))

; --- Heap declaration ---
     (insert "C ... Heap\n      integer w(1)\n      common /w/ w\n")

; ... Strip trailing blanks ...
     (goto-char (point-min)) (replace-regexp "[ \t]+$" "")

))

; --- macro fortran-find-cmdargs ---
(defun fortran-find-cmdargs (thisbuf varfilnam substart subend &optional arg)
"   Kernel to find command-line switches called by fortran-list-subroutine-args (&optional arg)"

   (interactive)

   (let (subnam cmdstr (pointnow (point))
	 (prefix-is-C-u
 		  (and arg (not (equal arg (prefix-numeric-value arg))))))

     ;; ... Handle each each occurence of a cmdopt(..)
;    (snit)
     (while (re-search-forward "cmdopt *('\\([a-zA-Z0-9_:.,+-]+=*\\)" subend t)
       (beginning-of-line)
       (if (looking-at "[Cc*]") (end-of-line)
	 (end-of-line)
	 (setq cmdstr (buffer-substring-no-properties (match-beginning 1) (match-end 1)))
	 (pop-to-buffer "*cmdlst*")
	 (insert (concat cmdstr " "))
	 (pop-to-buffer thisbuf)
	 )) ; end of while loop
     (goto-char pointnow)

     ;; ... Recursively search for command line arguments in subroutines called by this one
     (if (not arg) nil
;      (snit)
       (while (fortran-next-call subend) ; for each call, do
	 (setq subnam (buffer-substring-no-properties (match-beginning 1) (match-end 1)))
         ;; Handle cases where we don't follow tags: subnam -> nil
	 (if (not (is-tag subnam)) (setq subnam nil) ; skip if missing from tags table
	   (setq tagfile 			     ; file name containing tag
		 (save-excursion
		   (visit-tags-table-buffer) (file-of-tag)))
	   (if (string-match "slatsm" tagfile) (setq subnam nil)) ; skip files in slatsm directory
	   )
         ;; If subnam still not nil, tag should be acceptable; recursively find tags
	 (if (not subnam) nil
	   (setq pointnow (point))
	   (find-tag subnam)
	   (fortran-find-cmdargs (current-buffer) varfilnam (point) (save-excursion (end-of-fortran-subprogram) (point)) arg)
	   (pop-to-buffer thisbuf)
	   (goto-char pointnow)
	   )
	 ) ; end of while loop
       )

;     (snit)
))

; --- macro fortran-find-strux-elements ---
(defun fortran-find-strux-elements (thisbuf founds varfilnam strxnam substart subend &optional arg)
"   Kernel called by fortran-list-subroutine-args (&optional arg)"

   (interactive)

   (let (subnam cmdstr strxelt ifstore (pointnow (point))
	 (prefix-is-C-u
 		  (and arg (not (equal arg (prefix-numeric-value arg))))))

     ;; ... For each each occurence where an element is packed or unpacked, or passed as argument,
     ;;     add to declarations.
     (while (re-search-forward (concat "\\b" strxnam " *\\([(%]\\)" ) subend t)
       (if (save-excursion (beginning-of-line) (looking-at "[Cc*]")) nil
       (if (not (if (string= (buffer-substring-no-properties (match-beginning 1) (match-end 1)) "%")
	       (looking-at " *\\([A-za-z0-9_]+\\)")
	       (backward-char 1) (forward-sexp 1) (skip-chars-forward " ")
	       (if (not (looking-at "%")) nil (forward-char 1) (looking-at " *\\([A-za-z0-9_]+\\)"))))
	   nil ; nothing to do --- skip this reference

;	 (snit)


	 (setq strxelt (buffer-substring-no-properties (match-beginning 1) (match-end 1)))
	 (goto-char (match-end 1)) (if (looking-at " *(") (forward-sexp 1)) (skip-chars-forward " ")
					; character after element
	 ;; ifstore = "="  => element assigned a value
	 ;; ifstore = "[,)]"  => element passed as argument
	 (setq ifstore (buffer-substring-no-properties (point) (+ 1 (point))))


	 (pop-to-buffer "*arglst*")
					; First occurrence for this structure; set up comments
	 (if (not founds)
	     (progn (setq founds (point))
		    (end-of-line)
		    (insert "\n  Elts read: \n  Stored:    \n  Allocated: \n  Elts passed:\n  Passed to: ")))

					; Find 'Read' or 'Stored' apropos to I/O comments for this structure
	 (goto-char founds)
;	 (if (string-match "lshft" strxelt) (snit))
	 (if (string= ifstore "=")
	     (search-forward "Stored")
	   (if (or (string= ifstore ",") (string= ifstore ")"))
	       (search-forward "Elts passed")
	   (search-forward "Read")))
					; Append packlist
	 (end-of-line) (insert " " strxelt)
					; move past declaration of this structure
	 (beginning-of-line) (while (looking-at " ") (next-line 1))
	 (pop-to-buffer thisbuf)
	 ) ; end of insertion this strxelt
;     (snit)
       )) ; end of while loop

     ;; Calls to sitepack, spec2class, ptr_site
     (if (not (string= "s_site" strxnam)) nil
;       (snit)
       (goto-char substart)
       (while (re-search-forward "call ptr_site" subend t)
	 (if (save-excursion (beginning-of-line) (looking-at "[Cc*]")) nil
         (if (not (looking-at " *( *s_site *,[^,]+,'\\([A-za-z0-9_]+\\)'")) nil
	 (setq strxelt (buffer-substring-no-properties (match-beginning 1) (match-end 1)))
	 (pop-to-buffer "*arglst*")
	 (goto-char founds)
	 (search-forward "Allocated")
					; Append packlist
	 (end-of-line) (insert " " strxelt)
					; move past declaration of this structure
	 (beginning-of-line) (while (looking-at " ") (next-line 1))
	 (pop-to-buffer thisbuf))))

       (goto-char substart)
       (while (re-search-forward (concat "sitepack.+,'") subend t)
	 (if (save-excursion (beginning-of-line) (looking-at "[Cc*]")) nil
	 (if (looking-at "-") (progn (forward-char 1) (setq ifstore "=")) (setq ifstore nil))
	 (setq strxelt (buffer-substring-no-properties (point) (progn (skip-chars-forward "a-zA-Z0-9_") (point))))
	 (pop-to-buffer "*arglst*")
	 (goto-char founds)
	 (if (string= ifstore "=") (search-forward "Stored") (search-forward "Read"))
					; Append packlist
	 (end-of-line) (insert " " strxelt)
					; move past declaration of this structure
	 (beginning-of-line) (while (looking-at " ") (next-line 1))
	 (pop-to-buffer thisbuf)
	 )))

     (if (not (string= "s_spec" strxnam)) nil
       (goto-char substart)
       (while (re-search-forward "call ptr_spec" subend t)
	 (if (save-excursion (beginning-of-line) (looking-at "[Cc*]")) nil
         (if (not (looking-at " *( *s_spec *,[^,]+,'\\([A-za-z0-9_]+\\)'")) nil
	 (setq strxelt (buffer-substring-no-properties (match-beginning 1) (match-end 1)))
	 (pop-to-buffer "*arglst*")
	 (goto-char founds)
	 (search-forward "Allocated")
					; Append packlist
	 (end-of-line) (insert " " strxelt)
					; move past declaration of this structure
	 (beginning-of-line) (while (looking-at " ") (next-line 1))
	 (pop-to-buffer thisbuf))))

       (goto-char substart)
       (while (re-search-forward (concat "spec2class.+,'") subend t)
	 (if (save-excursion (beginning-of-line) (looking-at "[Cc*]")) nil
	 (if (looking-at "-") (progn (forward-char 1) (setq ifstore "=")) (setq ifstore nil))
	 (setq strxelt (buffer-substring-no-properties (point) (progn (skip-chars-forward "a-zA-Z0-9_") (point))))
	 (pop-to-buffer "*arglst*")
	 (goto-char founds)
	 (if (string= ifstore "=") (search-forward "Stored") (search-forward "Read"))
					; Append packlist
	 (end-of-line) (insert " " strxelt)
					; move past declaration of this structure
	 (beginning-of-line) (while (looking-at " ") (next-line 1))
	 (pop-to-buffer thisbuf)
	 )))

     (if (not (string= "s_array" strxnam)) nil
;       (snit)
       (goto-char substart)
       (while (re-search-forward "call ptr_arr" subend t)
	 (if (save-excursion (beginning-of-line) (looking-at "[Cc*]")) nil
         (if (not (looking-at " *( *s_array *,[^,]+,'\\([A-za-z0-9_]+\\)'")) nil
	   (setq strxelt (buffer-substring-no-properties (match-beginning 1) (match-end 1)))
	   (pop-to-buffer "*arglst*")
	   (goto-char founds)
	   (search-forward "Allocated")
					; Append packlist
	   (end-of-line) (insert " " strxelt)
					; move past declaration of this structure
	   (beginning-of-line) (while (looking-at " ") (next-line 1))
	   (pop-to-buffer thisbuf))))
       )

     (if (not (string= "s_bz" strxnam)) nil
;       (snit)
       (goto-char substart)
       (while (re-search-forward "call ptr_bz" subend t)
	 (if (save-excursion (beginning-of-line) (looking-at "[Cc*]")) nil
         (if (not (looking-at " *( *s_bz *,[^,]+,'\\([A-za-z0-9_]+\\)'")) nil
	   (setq strxelt (buffer-substring-no-properties (match-beginning 1) (match-end 1)))
	   (pop-to-buffer "*arglst*")
	   (goto-char founds)
	   (search-forward "Allocated")
					; Append packlist
	   (end-of-line) (insert " " strxelt)
					; move past declaration of this structure
	   (beginning-of-line) (while (looking-at " ") (next-line 1))
	   (pop-to-buffer thisbuf))))
       )

     (if (not (string= "s_lat" strxnam)) nil
;       (snit)
       (goto-char substart)
       (while (re-search-forward "call ptr_lat" subend t)
	 (if (save-excursion (beginning-of-line) (looking-at "[Cc*]")) nil
         (if (not (looking-at " *( *s_lat *,[^,]+,'\\([A-za-z0-9_]+\\)'")) nil
	   (setq strxelt (buffer-substring-no-properties (match-beginning 1) (match-end 1)))
	   (pop-to-buffer "*arglst*")
	   (goto-char founds)
	   (search-forward "Allocated")
					; Append packlist
	   (end-of-line) (insert " " strxelt)
					; move past declaration of this structure
	   (beginning-of-line) (while (looking-at " ") (next-line 1))
	   (pop-to-buffer thisbuf))))
       )

     (if (not (string= "s_ham" strxnam)) nil
;       (snit)
       (goto-char substart)
       (while (re-search-forward "call ptr_ham" subend t)
         (if (not (looking-at " *( *s_ham *,[^,]+,'\\([A-za-z0-9_]+\\)'")) nil
	   (setq strxelt (buffer-substring-no-properties (match-beginning 1) (match-end 1)))
	   (pop-to-buffer "*arglst*")
	   (goto-char founds)
	   (search-forward "Allocated")
					; Append packlist
	   (end-of-line) (insert " " strxelt)
					; move past declaration of this structure
	   (beginning-of-line) (while (looking-at " ") (next-line 1))
	   (pop-to-buffer thisbuf)))
       )

     (if (not (string= "s_pot" strxnam)) nil
;       (snit)
       (goto-char substart)
       (while (re-search-forward "call ptr_pot" subend t)
	 (if (save-excursion (beginning-of-line) (looking-at "[Cc*]")) nil
         (if (not (looking-at " *( *s_pot *,[^,]+,'\\([A-za-z0-9_]+\\)'")) nil
	   (setq strxelt (buffer-substring-no-properties (match-beginning 1) (match-end 1)))
	   (pop-to-buffer "*arglst*")
	   (goto-char founds)
	   (search-forward "Allocated")
					; Append packlist
	   (end-of-line) (insert " " strxelt)
					; move past declaration of this structure
	   (beginning-of-line) (while (looking-at " ") (next-line 1))
	   (pop-to-buffer thisbuf))))
       )

     (if (not (string= "s_str" strxnam)) nil
;       (snit)
       (goto-char substart)
       (while (re-search-forward "call ptr_str" subend t)
	 (if (save-excursion (beginning-of-line) (looking-at "[Cc*]")) nil
         (if (not (looking-at " *( *s_str *,[^,]+,'\\([A-za-z0-9_]+\\)'")) nil
	   (setq strxelt (buffer-substring-no-properties (match-beginning 1) (match-end 1)))
	   (pop-to-buffer "*arglst*")
	   (goto-char founds)
	   (search-forward "Allocated")
					; Append packlist
	   (end-of-line) (insert " " strxelt)
					; move past declaration of this structure
	   (beginning-of-line) (while (looking-at " ") (next-line 1))
	   (pop-to-buffer thisbuf))))
       )

     ;; Collect calls to subroutines that pass this structure as an argument
     (if (not founds) nil
;       (snit)
       (goto-char substart) (end-of-line)
       (while (re-search-forward (concat strxnam "\\([,)]\\)") subend t)
					; if at end of call, put point inside call arguments
	 (if (string= (buffer-substring-no-properties (match-beginning 1) (match-end 1)) ")") (backward-char 1))
					; skip past anything commented out
	 (if (save-excursion (beginning-of-line) (looking-at "[Cc*]")) nil
	   (save-excursion
	     (backward-up-list 1) (skip-chars-backward "[ \t]")
	     (setq subnam (buffer-substring-no-properties (point) (progn (skip-chars-backward "[A-Za-z0-9_]") (point))))
	     (backward-word 1)
	     (if (looking-at "\\(function\\|subroutine\\|SUBROUTINE\\)") (setq subnam nil))
;       (snit)
	     (if (and subnam (string-match "ptr_" subnam)) (setq subnam nil))
	     (if (string= subnam "sitepack") (setq subnam nil))
	     (if (string= subnam "spec2class") (setq subnam nil))

	     (beginning-of-line)
	     (if (looking-at "[Cc*]") (setq subnam nil)))
					; Poke list into *arglst*
	   ;; ... Recursively search for arguments in subroutines called by this one
	   (if (not subnam) nil
	     (pop-to-buffer "*arglst*")
	     (goto-char founds)
	     (search-forward "Passed to:") (end-of-line) (insert (concat " " subnam))

	     ;; If subnam still not nil, tag should be acceptable; recursively find tags
	     (if (not arg) nil
	       (pop-to-buffer thisbuf)
	       (if (not (is-tag subnam)) (setq subnam nil) ; skip if missing from tags table
		 (setq tagfile 			           ; file name containing tag
		       (save-excursion
			 (visit-tags-table-buffer) (file-of-tag)))
		 (if (string-match "slatsm" tagfile) (setq subnam nil)) ; skip files in slatsm directory
		 )
	       (if (not subnam) nil ;; do recursive serach
		 (setq pointnow (point))
		 (find-tag subnam)
;   	   (snit)
		 (fortran-find-strux-elements (current-buffer) founds varfilnam strxnam (point) (save-excursion (end-of-fortran-subprogram) (point)) arg)
;   	   (snit)
		 (pop-to-buffer thisbuf)
		 (goto-char pointnow)
		 ) ; end of recursive search
	     ))

	   (pop-to-buffer thisbuf)
	   ))
       ) ; (if (not founds) ...
     ) ; let
;    (snit)
)
