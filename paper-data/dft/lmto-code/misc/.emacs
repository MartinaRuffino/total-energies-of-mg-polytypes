; --- Miscellaneous setup ---
(display-time)
(setq abbrev-mode t)
(setq auto-save-default t)
(setq auto-save-interval 600)
(defvar fortran-comment-region "C")
;(defvar fortran-continuation-char ?.)

; --- Additional bindings ---
;(global-unset-key "\C- ")
;(global-set-key "\C-@"   'set-mark-command)
(global-set-key "\e\C-r" 'isearch-backward-regexp)
(global-set-key "\e."    'find-tagx)
(global-set-key "\e}"    'up-list)
(global-set-key ""     'isearch-forward)
(global-set-key ""   'save-buffer)
(global-set-key "\C-xv"  'insert-loop-var)
(global-set-key "\C-h\C-h" 'delete-backward-char)
(global-set-key "\C-cl"  'load-file)
(global-set-key "\C-xd"  'kill-word-surrounding-point)
;(global-set-key "S-down"  'next-comment-line)
(global-unset-key "\C-x\C-d")

; --- Mode-specific bindings ---
(switch-to-buffer "*scratch*") (fortran-mode)
(define-key fortran-mode-map "\C-c\c"   'fortran-comment-region)
(define-key fortran-mode-map "\C-co"    'fortran-var-as-offset)
(define-key fortran-mode-map "\C-cw"    'fortran-compare-progs)
(define-key fortran-mode-map "\C-c\C-o" 'fortran-offset-as-var)
(define-key fortran-mode-map "\C-c\C-f" 'fortran-forward-whitespace)
(define-key fortran-mode-map "\C-cf"    'fortran-forward-variable)
(define-key fortran-mode-map "\C-c\C-b" 'fortran-backward-whitespace)
(define-key fortran-mode-map "\C-cb"    'fortran-backward-variable)
(define-key fortran-mode-map "\en"      'fortran-next-block)
(define-key fortran-mode-map "\ep"      'fortran-previous-block)
(define-key fortran-mode-map "\C-c\\"   'fortran-collapse-continuation-line)
(define-key fortran-mode-map "\C-c\C-a" 'fortran-first-executable-statement)
(define-key fortran-mode-map "\C-cl"    'fortran-list-output-variables)
(define-key fortran-mode-map "\C-c\C-k" 'fortran-kill-word)
(define-key occur-mode-map "\C-cl"       'goto-line-occur-mode)

(setq c-mode-hook
      '(lambda ()
	(setq comment-line-regex
         "/\* *---")))

(setq fortran-mode-hook
      '(lambda ()
	(make-local-variable 'fortran-mode) (setq fortran-mode t)
	(setq comment-line-regex2 "^[Cc] *\\(--+\\)\\|\\(\\.\\.+\\.+\\)[ ()a-zA-Z0-9]\\|^ +subroutine")
	(setq comment-line-regex "^[Cc] *--+[ ()a-zA-Z0-9]\\|^ +subroutine")))

; --- Setup for TeX ---
(setq TeX-mode-hook
      '(lambda () (abbrev-mode 1)
	(setq paragraph-start "^[ \t]*$\\|^[ \t]*\\$\\$")
	(setq paragraph-separate "^[ \t]*$\\|\\$\\$[ \t]*$")
	(auto-fill-mode 1) (set-fill-column 65)
	(setq comment-line-regex "^% *---\\|^\\\\\\(hd\\)")
	(local-set-key "\C-ct"   'tex-macro-expand)
	(local-set-key "\C-cc"   'tex-comment-region)))
(setq LaTeX-mode-hook
      '(lambda ()
	(setq paragraph-start
	 "^[ \t]*$\\|^\\\\begin\\|^\\\\item\\|^[ \t]*\\$\\$")
	(setq paragraph-separate
	 "^[ \t]*$\\|^\\\\begin\\|^\\\\item\\|\\$\\$[ \t]*$")
	(setq comment-line-regex
	 "^% *---\\|^\\\\\\(section\\|subsection\\|begin\\)")
	(local-set-key "\C-ct"   'tex-macro-expand)))

(setq text-mode-hook
      '(lambda ()
	(setq comment-line-regex "^[#]* *---")))

(defvar tex-comment-region "%"
  "*String inserted by \\[fortran-comment-region] at start of each line in region.")

; --- macro tex-macro-expand ---
(defun tex-macro-expand ()
" Expand next Tex macro"
   (interactive)

   (while (re-search-forward "^\\\\def" (point-max) t)
     (save-excursion (progn
	(setq dnam
	      (buffer-substring 
		(point) (progn (forward-sexp 1) (point))))
	 (setq snam "")
	 (setq pm (save-excursion (progn (forward-sexp 1) (point))))
	 (while (< (point) pm)
	   (if (looking-at "[\\]")
	       (setq snam (progn (forward-char 1) (concat snam "\\\\")))
	       (setq snam (concat snam (buffer-substring (point)
				 (progn (forward-char 1) (point)))))))
	 (query-replace-regexp (concat "\\" dnam "\\b") snam)
       )))
   (query-replace-regexp "}\\\\ "  "} ")
   (query-replace-regexp "}\\\\$" "}" nil)
)

; --- macro tex-comment-region ---
(defun tex-comment-region (beg-region end-region arg)
  "Comments every line in the region.
Puts tex-comment-region at the beginning of every line in the region.
BEG-REGION and END-REGION are args which specify the region boundaries.
With non-nil ARG, uncomments the region."
  (interactive "*r\nP")
  (let ((end-region-mark (make-marker)) (save-point (point-marker)))
    (set-marker end-region-mark end-region)
    (goto-char beg-region)
    (beginning-of-line)
    (if (not arg)			;comment the region
	(progn (insert tex-comment-region)
	       (while (and  (= (forward-line 1) 0)
			    (< (point) end-region-mark))
		 (insert tex-comment-region)))
      (let ((com (regexp-quote tex-comment-region))) ;uncomment the region
	(if (looking-at com)
	    (delete-region (point) (match-end 0)))
	(while (and  (= (forward-line 1) 0)
		     (< (point) end-region-mark))
	  (if (looking-at com)
	      (delete-region (point) (match-end 0))))))
    (goto-char save-point)
    (set-marker end-region-mark nil)
    (set-marker save-point nil)))

; --- Create shell ---
;(load-file "~/shell.el")
(setq shell-pushd-regexp "p\\(ush\\)*d")
(shell) (local-set-key "c" (quote clear-buffer))
(switch-to-buffer "*scratch*") (kill-buffer "*scratch*")
;(define-key shell-mode-map "\C-cr"       'comint-previous-matching-input-from-input)

; --- Keyboard macros ---
(fset 'get-semi-bulk-mod
   "dval '4/9*bbfb\\z /bb>\\*147'")

(fset 'get-shear-mod
   "<  k k")

(fset 'apollo-kd-DM-buffer
   "b*paste*kb*paste*kd f9 es 'HE' ke HDk")

(fset 'copy-region-from-x
   "xinsert-registerx")

(fset 'copy-region-to-x
   "xcopy-to-registerx")

(fset 'set-word-for-isearch
   "[a-z9-0]*xisearch-fo	-re	\\b\\b")

(fset 'goto-line-occur-mode
   "10")

(fset 'clear-buffer
   "xmark-whole-bufferxkill-region")

(fset 'find-function-args
   ".( wb")
(global-set-key "\C-ca" 'find-function-args)

(fset 'fortran-collapse-continuation-line
   "\\")

(fset 'recenter-this-line-at-top
   "\C-u0\C-l")

;(fset 'sort-lmovl
;   "ovl lmo10xsort-fields lmo")

(fset 'test
   "l~/t.el
xdebug
f-x
")
(global-set-key "\C-c\C-t" 'test)
(global-set-key "\C-ct" 'fortran-call-tree)

; --- macro line-len ---
(defun line-len ()
"  Returns length of line point is on."
   (interactive)
   (- (save-excursion (end-of-line) (point)) (save-excursion (beginning-of-line) (point))))

(defvar comment-line-regex "\\(^[Cc;%#]* *---\\)\\|\\(---[ ']*$\\)"
  "regular expression for comment-line macros")
;(defvar comment-line-regex "^[Cc;%]* *---"
;  "regular expression for comment-line macros")
(make-variable-buffer-local 'comment-line-regex)

; --- macro list-comment-lines (25 Nov 94) ---
(defun list-comment-lines (&optional arg)
" lists in a window lines matching variable 'comment-line-regex.'
 Optional ARG does the following:
  ARG = C-U use 'comment-line-regex2' instead of 'comment-line-regex'
  ARG = C-U C-U prompts for string to use for search."

  (interactive "P")

  (let (strn
	(prefix-is-C-u
	 (and arg (not (equal current-prefix-arg
			  (prefix-numeric-value current-prefix-arg))))))

    (if prefix-is-C-u
;... Non-numeric ARG
	(if (= 16 (prefix-numeric-value current-prefix-arg))
	  (setq strn (read-from-minibuffer "match string: "))
	  (setq strn comment-line-regex2))
;... Numeric or no ARG
      (setq strn comment-line-regex))

    (save-excursion (beginning-of-buffer) (list-matching-lines strn))
    (pop-to-buffer "*Occur*")
;   (delete-other-windows-quietly)
))

; --- macro previous-comment-line (25 Nov 94) ---
(defun previous-comment-line (&optional arg)
" Moves up to previous comment line, putting line at top of window.
 Comment line defined by variable 'comment-line-regex.'
 Optional ARG does the following:
  ARG = C-U use 'comment-line-regex2' instead of 'comment-line-regex'
  ARG = C-U C-U prompts for string to use for search
  ARG = anything else, moves up ARG lines."

  (interactive "P")
  (let (strn
	(prefix-is-C-u
	 (and arg (not (equal current-prefix-arg
			  (prefix-numeric-value current-prefix-arg))))))

    (if prefix-is-C-u
;... Non-numeric ARG
      (if (= 16 (prefix-numeric-value current-prefix-arg))
	  (setq strn (read-from-minibuffer "search string: ") nsrch 1)
	  (setq strn comment-line-regex2) nsrch 1)
;... Numeric or no ARG
      (progn
	(setq strn comment-line-regex nsrch 1)
	(if arg (setq strn comment-line-regex2
		      nsrch (prefix-numeric-value arg)))))

    (if (= (point) (point-min)) (progn (ding) (goto-char (point-max))))

    (recenter 0)
    (if (re-search-backward strn nil t nsrch)
	(progn (beginning-of-line 1) (recenter 0))
	(progn (goto-char (point-min)) (recenter 0))))
)

; --- macro next-comment-line (25 Nov 94) ---
(defun next-comment-line (&optional arg)
" Moves down to next comment line, putting line at top of window.
 Comment line defined by variable 'comment-line-regex.'
 Optional ARG does the following:
  ARG = C-U use 'comment-line-regex2' instead of 'comment-line-regex'
  ARG = C-U C-U prompts for string to use for search
  ARG = anything else, moves down ARG lines."

  (interactive "P")
  (let (strn
	(prefix-is-C-u
	 (and arg (not (equal current-prefix-arg
			  (prefix-numeric-value current-prefix-arg))))))

    (if prefix-is-C-u
;... Non-numeric ARG
      (if (= 16 (prefix-numeric-value current-prefix-arg))
	  (setq strn (read-from-minibuffer "search string: ") nsrch 1)
	  (setq strn comment-line-regex2) nsrch 1)
;... Numeric or no ARG
      (progn
	(setq strn comment-line-regex nsrch 1)
	(if arg (setq strn comment-line-regex2
		      nsrch (prefix-numeric-value arg)))))

    (if (= (point) (point-max)) (progn (ding) (goto-char (point-min))))

    (recenter 0)
    (move-to-window-line 1)
    (if (re-search-forward strn nil t nsrch)
	(progn (beginning-of-line 1) (recenter 0))
	(progn (move-to-window-line 0) (goto-char (point-max)))))
)

; --- macro kill-word-surrounding-point ---
(defun kill-word-surrounding-point (&optional arg)
" Like kill-word, but kills entire word surrounding point when inside.
  a word.  With argument, kill arg times."
  (interactive "p") (if (not arg) (setq arg 1))

  (if (looking-at "\\w") (progn (forward-word 1) (backward-word 1)))
  (kill-word arg)
)
; --- macro what-cursor-position ---
(defun what-cursor-position ()
  "Print info on cursor position (on screen and within buffer).
   This version modified to display line number and max column"
  (interactive)
  (let* ((char (following-char))
	 (beg (point-min))
	 (end (point-max))
         (pos (point))
	 (total (buffer-size))
         (eol (line-len))
	 (cln (count-lines (point-min) (point)))
	 (mxln (count-lines 1 (point-max)))
	 (percent (if (> total 50000)
		      ;; Avoid overflow from multiplying by 100!
		      (/ (+ (/ total 200) (1- pos)) (max (/ total 100) 1))
		    (/ (+ (/ total 2) (* 100 (1- pos))) (max total 1))))
	 (hscroll (if (= (window-hscroll) 0)
		      ""
		    (format " Hscroll=%d" (window-hscroll))))
	 (col (current-column)))
    (if (= (current-column) 0) (setq cln (+ cln 1)))
    (if (= pos end)
	(if (or (/= beg 1) (/= end (1+ total)))
	    (message "point=%d of %d(%d%%) <%d - %d>  col %d(%d) line %d(%d) %s"
		     pos total percent beg end col eol cln mxln hscroll)
	  (message "point=%d of %d(%d%%)   col %d(%d) line %d(%d) %s"
		   pos total percent col eol cln mxln hscroll))
      (if (or (/= beg 1) (/= end (1+ total)))
	  (message "Char: %s (0%o)  point=%d of %d(%d%%) <%d - %d>  column %d %s"
		   (single-key-description char) char pos total percent beg end col hscroll)
	(message "Char: %s (0%o)  point=%d of %d(%d%%)  col %d(%d) line %d(%d) %s"
		 (single-key-description char) char pos total percent col eol cln mxln hscroll)))))

; --- macro count-words-region ---
(defun count-words-region (start end)
  "Print number of lines,words,chars in the region."
  (interactive "r")
  (message "Region has %d lines, %d words, %d characters"
   (count-lines start end) (count-words start end) (- end start)))

(defun count-words (start end)
  "Return number of words between START and END."
  (save-excursion
    (save-restriction
      (narrow-to-region start end)
      (goto-char (point-min)) (setq count 0)
      (while (forward-word 1) (setq count (+ count 1)))))
  count)

; --- macro fortran-offset-as-var ---
;(fset 'fortran-offset-as-var
;   "xfortran-forward-va	xset-mark-co	xfort	backward-va	")
(fset 'fortran-offset-as-var
   "xup-listxset-mark-com	-1xset-mark-com	")

; --- macro fortran-kill-word ---
(defun fortran-kill-word (&optional arg)
"Kill the next fortran word.
With argument, do this that many times."
  (interactive "p") (if (not arg) (setq arg 1))

  (kill-region
    (point) (save-excursion (fortran-forward-variable arg) (point)))

)

; --- macro fortran-var-as-offset ---
(defun fortran-var-as-offset ()
"  Converts variable zz to w(ozz)"
   (interactive)

   (if (not (looking-at "\\<")) (fortran-backward-variable))
   (insert-string "w(o") (fortran-forward-variable) (insert-string ")")
)

; --- macro compare-buffers ---
(defun compare-buffers (&optional arg)
" Find first character common to two windows with minimum change
  in buffer position, and execute compare-windows.
  ding sounds if case mismatch, but compare-windows continues.
  Optional ARG moves point in each window before comparing windows:
    ARG > 0   Moves forward ARG characters
    ARG = C-U Moves point past next whitespace
    ARG = C-U C-U Moves point to beginning of next line
    ARG = C-U C-U C-U Prompts for regex, Moves point past regex
    ARG < 0   Moves point down -ARG lines"

  (interactive "P")

  (let (p1 p2 b1 b2 w1 w2 charmatch
	totskip skip1 skip2 skip1x skip2x char1 char2 regex
	(prefix-is-C-u
	 (and arg (not (equal current-prefix-arg
			 (prefix-numeric-value current-prefix-arg))))))

;   ... arg = 0 if no prefix, -1 if C-u, otherwise numerical value
    (if prefix-is-C-u
	(if (= 16 (prefix-numeric-value current-prefix-arg))
	    (setq arg -1)
	    (setq arg 0))
	(if arg (setq arg (prefix-numeric-value arg))
	    (setq arg 0)))

    (setq w1 (selected-window))
    (setq w2 (next-window (selected-window)))
    (if (eq w2 (selected-window)) (error "No other window."))

;   ... Move forward ARG characters or -ARG lines if ARG < 0
    (if (>= arg 0)
	(if (> arg 0)
;           ... Case ARG>0: move forward ARG characters
	    (progn (goto-char (+ (point) arg))
		   (set-window-point w2 (+ (window-point w2) arg)))
;           ... Case ARG = C-U: move past next whitespace
	    (if (= 4 (prefix-numeric-value current-prefix-arg))
		(progn
		 (select-window w2) (re-search-forward "[\t\n ]+" nil 0)
		 (set-window-point w2 (point))
		 (select-window w1) (re-search-forward "[\t\n ]+" nil 0)
		 (set-window-point w1 (point))))
;           ... Case ARG = C-U C-U C-U: search both windows regex
	    (if (= 64 (prefix-numeric-value current-prefix-arg))
		(progn
		  (setq regex (read-string "RE search: "))
		 (select-window w2) (re-search-forward regex nil 0)
		 (set-window-point w2 (point))
		 (select-window w1) (re-search-forward regex nil 0)
		 (set-window-point w1 (point))))
	    )
;       ... Case ARG<0: move down -ARG lines
	(progn
	  (select-window w2) (beginning-of-line) (previous-line arg)
	  (set-window-point w2 (point))
	  (select-window w1) (beginning-of-line) (previous-line arg)
	  (set-window-point w1 (point))))

;   ... Not doing this messes up emacs' real p2, in window 2
    (other-window 1) (other-window -1)
    (setq p1 (point) b1 (current-buffer))
    (setq p2 (window-point w2)  b2 (window-buffer w2))

;   ... Setup for loop to force char match
    (setq char1 (char-after (point)))
    (set-buffer b2) (setq char2 (char-after (window-point w2)))
    (setq totskip 0)
;   ... In this loop, either totskip increments or char1=char2
    (while (not (char-equal char1 char2))
      (setq totskip (+ totskip 1))
      (setq skip1 0)
      (while (<= skip1 totskip)
	(set-buffer b1)
	(setq charmatch (char-after (+ skip1 (point))))
	(setq skip2 0)
	(while (<= (+ skip1 skip2) totskip)
	  (set-buffer b2)
	  (if (char-equal charmatch (char-after (+ skip2 (point))))
	      (progn
		(setq char1 charmatch char2 charmatch
		      skip1x skip1 skip2x skip2)
		(set-window-point w2 (+ p2 skip2))
		(set-buffer b1) (goto-char (+ p1 skip1))
		(setq skip1 totskip)))
	  (setq skip2 (+ skip2 1)))
	(setq skip1 (+ skip1 1))))

    (set-buffer b1)
    (if (not (= totskip 0))
	(message "skipping %d,%d chars ..." skip1x skip2x))

;   ... Reset char1, char2
    (setq char1 (char-after (point)))
    (save-excursion
      (set-buffer b2) (setq char2 (char-after (window-point w2))))

;   ... compare-windows, allowing for different case
    (while (char-equal char1 char2)
      (compare-windows)
      (setq char1 (char-after (point)))
      (save-excursion
	(set-buffer b2)
	(setq char2 (char-after (window-point w2))))
      (if (and (char-equal char1 char2) (not (= char1 char2)))
	  (progn (goto-char (+ (point) 1))
		 (set-window-point w2 (+ (window-point w2) 1))
		 ))
      )
;   ... Not doing this messes up emacs' record of p2, in window 2
    (other-window 1) (other-window -1)
    )
)

; --- macro fortran-compare-progs ---
(defun fortran-compare-progs (&optional arg)
" Similar to compare-buffers, except tailored to fortran programs.
  fortran-forward-whitespace(0) is run in both buffers before
  compare-buffers(0) is called; optional ARG makes additional
  movements:
    ARG > 0   executes fortran-forward-variable (ARG)
    ARG = C-U executes fortran-next-statement
    ARG < 0   executes fortran-next-statement -ARG times
"

  (interactive "P")

  (let
      ((prefix-is-C-u
	(and arg (not (equal current-prefix-arg
		         (prefix-numeric-value current-prefix-arg)))))
       w1 w2)

; arg = 0 if no prefix, -1 if C-u, otherwise numerical value
    (if prefix-is-C-u
	(setq arg -1)
	(if arg (setq arg (prefix-numeric-value arg))
	    (setq arg 0)))

    (setq w1 (selected-window))
    (setq w2 (next-window (selected-window)))

; Run fortran-forward-whitespace
    (select-window w2) (fortran-forward-whitespace 0)
    (select-window w1) (fortran-forward-whitespace 0)

; Move forward ARG characters or -ARG lines if ARG < 0
  (if (>= arg 0)
      (if (> arg 0)
	  (progn
	    (select-window w2) (fortran-forward-variable arg)
	    (select-window w1) (fortran-forward-variable arg)))
	  (while (< arg 0)
	    (select-window w2) (fortran-next-statement)
	    (select-window w1) (fortran-next-statement)
	    (setq arg (+ arg 1))))

; compare-buffers
  (compare-buffers 0)

))

; --- macro select-last-buffer ---
(defun select-last-buffer (&optional arg)
" Selects last buffer from list in other window.
 Optional ARG moves selects ARGth-to-last buffer"
  (interactive "p") (if (not arg) (setq arg 1))

  (let ((list (buffer-list)))
    (while (not (= arg 0)) (setq arg (- arg 1)) (setq list (cdr list)))
    (setq buffer (car list)) (setq bnam (buffer-name buffer))
    (if (string= bnam " *Minibuf-0*") (setq list (cdr list)))
    (setq buffer (car list)) (setq bnam (buffer-name buffer))
    (switch-to-buffer buffer))
)

; --- macro select-last-buffer-1-window ---
(defun select-last-buffer-1-window (&optional arg)
" Selects last buffer from list in one window.
 Optional ARG moves selects ARGth-to-last buffer"

  (interactive "P")
  (if (not arg) (setq arg 1))

  (select-last-buffer arg)
  (delete-other-windows-quietly)
)

(defvar *loop-var-now 0
  "is inserted as a string when 'insert-loop-var is invoked.")
(defvar *loop-var-inc 1
  "Amount *loop-var-now is incremented when 'insert-loop-var is invoked.")
(defvar *loop-var-pad 0
 "*Minimum number of characters inserted when 'insert-loop-var invoked.")

; --- macro insert-loop-var ---
(defun insert-loop-var (&optional arg)
"Inserts variable '*loop-var-now' into current buffer and increments
*loop-var-now by *loop-var-inc.  Optional ARG does the following:
  ARG = C-U sets *loop-var-now = 0 and *loop-var-inc = 1
       (inserts nothing into current buffer)
  ARG = C-U C-U does same, but prompts for *loop-var and *loop-var-inc,
        and also variable *loop-var-pad. (Setting *loop-var-pad nonzero
        ensures that when *loop-var-now is inserted into the
        current buffer, at least *loop-var-pad characters are inserted.
        It does this by prepending leading zero's.)
  ARG = anything else, adds ARG - *loop-var-inc to *loop-var-now
        before inserting into buffer."

  (interactive "P")
  (let (strn
	(prefix-is-C-u
	 (and arg (not (equal current-prefix-arg
			  (prefix-numeric-value current-prefix-arg))))))

    (if prefix-is-C-u
;... Non-numeric ARG
	(if (= 16 (prefix-numeric-value current-prefix-arg))
	    (progn
	      (setq strn (read-string "Reset loop counter to: "))
	      (if (not (string= "" strn))
		  (setq *loop-var-now (string-to-int strn)))
	      (message "%s  %d" strn *loop-var-now)
	      (setq strn (read-string "Reset counter increment to: "))
	      (if (not (string= "" strn))
		  (setq *loop-var-inc (string-to-int strn)))
	      (message "%s  %d" strn *loop-var-inc)
	      (setq strn (read-string "Reset counter pad to: "))
	      (if (not (string= "" strn))
		  (setq *loop-var-pad (string-to-int strn))))
	    (progn
	      (setq *loop-var-now 0)
	      (setq *loop-var-inc 1)))
;... Numeric or no ARG
	(progn
	  (if arg (setq *loop-var-now
			(+ *loop-var-now (- arg *loop-var-inc))))
	  (setq strn (int-to-string *loop-var-now))
	  (if (> *loop-var-pad 0)
	      (while (< (length strn) *loop-var-pad)
		(setq strn (concat "0" strn))))
	  (insert strn)
	  (setq *loop-var-now (+ *loop-var-now *loop-var-inc))
))))

; --- macro clear-tabs-in-buffer ---
(defun clear-tabs-in-buffer ()
  (interactive) (untabify (point-min) (point-max)))

; --- macro split-window-quietly ---
(defun split-window-quietly (&optional arg)
  "Split the window vertically with minimum redisplay.
This window becomes the uppermost of the two, and gets
ARG lines.  No arg means split equally."
  (interactive "P")
  (let* (
	  (num-arg (and arg (prefix-numeric-value arg)))
	  (oldpt (point))
	  (oldstart (window-start))
	  (scroll (or num-arg (/ (window-height) 2)))
	  (scrollpoint (progn (move-to-window-line scroll)
			 (point)))
	  (barstart (progn
		      (move-to-window-line (- scroll 1))
		      (point)))
	  )
    (split-window nil num-arg)
    (goto-char oldstart)
    (recenter 0)
    (other-window 1)
    (goto-char scrollpoint)
    (recenter 0)
    (if (< oldpt scrollpoint )
      (if (>= oldpt barstart)
	(progn
	  (other-window -1)
	  (move-to-window-line (- scroll 2))
	  )
	(progn
	  (other-window -1)
	  (goto-char oldpt)
	  ))
      (progn
	(goto-char oldpt)
	))
    ))

; --- macro delete-other-windows-quietly ---
(defun delete-other-windows-quietly ()
  "Delete other windows with minimum redisplay"
  (interactive)
  (let* ((oldpt (point))
	  (oldtopchar (window-start))
	  (oldtop (car (cdr (window-edges))))
	  )
    (delete-other-windows (selected-window))
    (goto-char oldtopchar)
    (recenter oldtop)
    (goto-char oldpt)))

; --- macro backward-down-list ---
(defun backward-down-list (&optional arg)
"Move forward down one level of parentheses.
With argument, do this that many times.
A negative argument means move backward but still go down a level."
  (interactive "p") (if (not arg) (setq arg 1))
  (down-list (- 0 arg))
)

; --- macros remember-point, recover-point ---
(defun remember-point () (interactive) (point-to-register ?p))
(defun recover-point  () (interactive) (register-to-point ?p))

; -------- Function definitions for five special Buffers ------
(setq buffer1 "-")
(setq buffer2 "-")
(setq buffer3 "-")
(setq buffer4 "-")
(setq buffer5 "-")

(defun define-buffer-number (num)
" Assigns current buffer to special set 1-5."
  (interactive "Nassign buffer to number: ")
  (if (= num 1) (setq buffer1 (buffer-name)))
  (if (= num 2) (setq buffer2 (buffer-name)))
  (if (= num 3) (setq buffer3 (buffer-name)))
  (if (= num 4) (setq buffer4 (buffer-name)))
  (if (= num 5) (setq buffer5 (buffer-name)))
  (show-buffer-definitions) )

(defun goto-buffer1 () (interactive) (switch-to-buffer buffer1))
(defun goto-buffer2 () (interactive) (switch-to-buffer buffer2))
(defun goto-buffer3 () (interactive) (switch-to-buffer buffer3))
(defun goto-buffer4 () (interactive) (switch-to-buffer buffer4))
(defun goto-buffer5 () (interactive) (switch-to-buffer buffer5))

(defun show-buffer-file-name () (interactive)
    (message "%s" (buffer-file-name)))

(defun show-buffer-definitions () (interactive)
  (message "%s  %s  %s  %s  %s" buffer1 buffer2 buffer3 buffer4 buffer5 ))

; --- Customization of Fortran variables (fortran.el loaded) ---
(load-file "~/fortran.el")
(load-file "~/tags.el")
 (set-variable (quote fortran-continuation-indent) 2)
 (set-variable (quote fortran-if-indent) 2)
 (set-variable (quote fortran-do-indent) 2)
 (set-variable (quote fortran-line-number-indent) 5)
 (set-variable (quote fortran-continuation-char) 46)
 (set-variable (quote fortran-comment-region) "C")
 (set-variable (quote fortran-check-all-num-for-matching-do) t)

; --- GOLD key on ^Z ---
  (global-unset-key "")
  (define-key global-map "" 'GOLD-prefix) ; GOLD-prefix on ^z
  (setq GOLD-map (make-keymap))
  (fset 'GOLD-prefix GOLD-map)
  (defvar GOLD-map nil "Map for GOLD prefix on ^Z")

; --- Miscellaneous GOLD-key bindings ---
  (define-key GOLD-map "\e2"    'split-window-quietly)
  (define-key GOLD-map "\e1"    'delete-other-windows-quietly)
  (define-key GOLD-map "\en"    'picture-move-down)
  (define-key GOLD-map "\ep"    'picture-move-up)
  (define-key GOLD-map "%"    'query-replace-regexp)
  (define-key GOLD-map ","    'remember-point)
  (define-key GOLD-map "."    'recover-point)
  (define-key GOLD-map "["    'backward-up-list)
  (define-key GOLD-map "]"    'up-list)
  (define-key GOLD-map "{"    'backward-down-list)
  (define-key GOLD-map "}"    'down-list)
  (define-key GOLD-map "b"    'select-last-buffer-1-window)
  (define-key GOLD-map "c"    'fortran-compare-call-args)
  (define-key GOLD-map "f"    'apollo-find-file)
  (define-key GOLD-map "g"    'goto-line)
  (define-key GOLD-map "i"    'set-word-for-isearch)
  (define-key GOLD-map "-"    'recenter-this-line-at-top)
  (define-key GOLD-map "k"    'kill-rectangle)
  (define-key GOLD-map "l"    'list-comment-lines)
  (define-key GOLD-map "p"    'picture-mode)
  (define-key GOLD-map "s"    'shell)
  (define-key GOLD-map "t"    'clear-tabs-in-buffer)
  (define-key GOLD-map "w"    'compare-buffers)
  (define-key GOLD-map "x"    'copy-region-from-x)
  (define-key GOLD-map "y"    'yank-rectangle)
  (define-key GOLD-map "="    'count-words-region)
  (define-key GOLD-map "\C-b" 'bury-buffer)
  (define-key GOLD-map "\C-p" 'previous-comment-line)
  (define-key GOLD-map "\C-n" 'next-comment-line)
  (define-key GOLD-map "\C-x" 'copy-region-to-x)
  (define-key GOLD-map "\C-z" 'suspend-emacs)
  (define-key GOLD-map "\t"   'picture-tab-search)
  (define-key GOLD-map "\ep"  'apollo-kd-DM-buffer)
; (define-key GOLD-map "0"    'aix-add-0)
; ... Keymaps for 5 special buffers
  (define-key GOLD-map "d" 'define-buffer-number)
  (define-key GOLD-map "q" 'show-buffer-definitions) ; query defs
  (define-key GOLD-map "1" 'goto-buffer1) ; b1
  (define-key GOLD-map "2" 'goto-buffer2) ; b2
  (define-key GOLD-map "3" 'goto-buffer3) ; b3
  (define-key GOLD-map "4" 'goto-buffer4) ; b4
  (define-key GOLD-map "5" 'goto-buffer5) ; b5

; --- VT100 Key bindings ---
(defun vt100-key-bindings ()
" Bindings specific to VT100 terminal"
  (interactive)
; ... Special to VT100 simulation
  (global-set-key " "   (quote set-mark-command))
; (load-file "~/.emacs-100")
; (load-file "~/fortran.el")
)

; --- For the AIX ---
(defun aix-add-0 () (interactive)
  (replace-regexp " \\(-?\\)\\.\\([0-9]\\)" "\\10.\\2" nil))

; --- For the Apollo ---
(defvar *preempt-display-manager-bindings* t)

; --- Terminal-dependent setup ---
(if (eq window-system 'apollo)
    (progn (load-file "~/.emacs-apollo") (my-apollo-key-bindings)))
(if (string= (getenv "TERM") "vt100") (vt100-key-bindings))

(put 'eval-expression 'disabled nil)
