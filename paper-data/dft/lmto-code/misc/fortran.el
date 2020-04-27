;;; fortran.el --- Fortran mode for GNU Emacs

;; Copyright (c) 1986, 1993, 1994, 1995, 1997, 1998 Free Software Foundation, Inc.

;; Author: Michael D. Prange <prange@erl.mit.edu>
;; Maintainer: Dave Love <fx@gnu.org>
;; Keywords: languages

;; This file is part of GNU Emacs.

;; GNU Emacs is free software; you can redistribute it and/or modify
;; it under the terms of the GNU General Public License as published by
;; the Free Software Foundation; either version 2, or (at your option)
;; any later version.

;; GNU Emacs is distributed in the hope that it will be useful,
;; but WITHOUT ANY WARRANTY; without even the implied warranty of
;; MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;; GNU General Public License for more details.

;; You should have received a copy of the GNU General Public License
;; along with GNU Emacs; see the file COPYING.  If not, write to the
;; Free Software Foundation, Inc., 59 Temple Place - Suite 330,
;; Boston, MA 02111-1307, USA.

;;; Commentary:

;; This mode is documented in the Emacs manual.
;;
;; Note that it is for editing Fortran77 or Fortran90 fixed source
;; form.  For editing Fortran 90 free format source, use `f90-mode'
;; (f90.el).

;;; History:

;; Fortran mode was upgraded by Stephen A. Wood (saw@cebaf.gov).

;; We acknowledge many contributions and valuable suggestions by
;; Lawrence R. Dodd, Ralf Fassel, Ralph Finch, Stephen Gildea,
;; Dr. Anil Gokhale, Ulrich Mueller, Mark Neale, Eric Prestemon,
;; Gary Sabot and Richard Stallman.

;;; Code:

;; Todo:

;; * Tidy it all up!  (including renaming non-`fortran' prefixed
;;   functions).
;; * Implement insertion and removal of statement continuations in
;;   mixed f77/f90 style, with the first `&' past column 72 and the
;;   second in column 6.
;; * Support any other extensions to f77 grokked by GNU Fortran.
;; * Change fontification to use font-lock-syntactic-keywords for
;;   fixed-form comments.  (Done, but doesn't work properly with
;;   lazy-lock in pre-20.4.)

(require 'easymenu)

(defgroup fortran nil
  "Fortran mode for Emacs"
  :link '(custom-manual "(emacs)Fortran")
  :group 'languages)

(defgroup fortran-indent nil
  "Indentation variables in Fortran mode"
  :prefix "fortran-"
  :group 'fortran)

(defgroup fortran-comment nil
  "Comment-handling variables in Fortran mode"
  :prefix "fortran-"
  :group 'fortran)


;;;###autoload
(defcustom fortran-tab-mode-default nil
  "*Default tabbing/carriage control style for empty files in Fortran mode.
A value of t specifies tab-digit style of continuation control.
A value of nil specifies that continuation lines are marked
with a character in column 6."
  :type 'boolean
  :group 'fortran-indent)

;; Buffer local, used to display mode line.
(defcustom fortran-tab-mode-string nil
  "String to appear in mode line when TAB format mode is on."
  :type '(choice (const nil) string)
  :group 'fortran-indent)
(make-variable-buffer-local 'fortran-tab-mode-string)

(defcustom fortran-do-indent 3
  "*Extra indentation applied to DO blocks."
  :type 'integer
  :group 'fortran-indent)

(defcustom fortran-if-indent 3
  "*Extra indentation applied to IF blocks."
  :type 'integer
  :group 'fortran-indent)

(defcustom fortran-structure-indent 3
  "*Extra indentation applied to STRUCTURE, UNION, MAP and INTERFACE blocks."
  :type 'integer
  :group 'fortran-indent)

(defcustom fortran-continuation-indent 5
  "*Extra indentation applied to Fortran continuation lines."
  :type 'integer
  :group 'fortran-indent)

(defcustom fortran-comment-indent-style 'fixed
  "*How to indent comments.
nil forces comment lines not to be touched,
'fixed makes fixed comment indentation to `fortran-comment-line-extra-indent'
columns beyond `fortran-minimum-statement-indent-fixed' (for
`indent-tabs-mode' of nil) or `fortran-minimum-statement-indent-tab' (for
`indent-tabs-mode' of t), and 'relative indents to current
Fortran indentation plus `fortran-comment-line-extra-indent'."
  :type '(radio (const :tag "Untouched" nil) (const fixed) (const relative))
  :group 'fortran-indent)

(defcustom fortran-comment-line-extra-indent 0
  "*Amount of extra indentation for text within full-line comments."
  :type 'integer
  :group 'fortran-indent
  :group 'fortran-comment)

(defcustom comment-line-start nil
  "*Delimiter inserted to start new full-line comment."
  :type '(choice string (const nil))
  :group 'fortran-comment)

(defcustom comment-line-start-skip nil
  "*Regexp to match the start of a full-line comment."
  :type '(choice string (const nil))
  :group 'fortran-comment)

(defcustom fortran-minimum-statement-indent-fixed 6
  "*Minimum statement indentation for fixed format continuation style."
  :type 'integer
  :group 'fortran-indent)

(defcustom fortran-minimum-statement-indent-tab (max tab-width 6)
  "*Minimum statement indentation for TAB format continuation style."
  :type 'integer
  :group 'fortran-indent)

;; Note that this is documented in the v18 manuals as being a string
;; of length one rather than a single character.
;; The code in this file accepts either format for compatibility.
(defcustom fortran-comment-indent-char " "
  "*Single-character string inserted for Fortran comment indentation.
Normally a space."
  :type 'string
  :group 'fortran-comment)

(defcustom fortran-line-number-indent 1
  "*Maximum indentation for Fortran line numbers.
5 means right-justify them within their five-column field."
  :type 'integer
  :group 'fortran-indent)

(defcustom fortran-check-all-num-for-matching-do nil
  "*Non-nil causes all numbered lines to be treated as possible DO loop ends."
  :type 'boolean
  :group 'fortran)

(defcustom fortran-blink-matching-if nil
  "*Non-nil causes \\[fortran-indent-line] on ENDIF statement to blink on matching IF.
Also, from an ENDDO statement blink on matching DO [WHILE] statement."
  :type 'boolean
  :group 'fortran)

(defcustom fortran-continuation-string "$"
  "*Single-character string used for Fortran continuation lines.
In fixed format continuation style, this character is inserted in
column 6 by \\[fortran-split-line] to begin a continuation line.
Also, if \\[fortran-indent-line] finds this at the beginning of a line, it will
convert the line into a continuation line of the appropriate style.
Normally $."
  :type 'string
  :group 'fortran)

(defcustom fortran-comment-region "c$$$"
  "*String inserted by \\[fortran-comment-region] at start of each \
line in region."
  :type 'string
  :group 'fortran-comment)

(defcustom fortran-electric-line-number t
  "*Non-nil causes line number digits to be moved to the correct \
column as typed."
  :type 'boolean
  :group 'fortran)

(defvar fortran-column-ruler-fixed
  "0   4 6  10        20        30        40        5\
0        60        70\n\
\[   ]|{   |    |    |    |    |    |    |    |    \
\|    |    |    |    |}\n"
  "String displayed above current line by \\[fortran-column-ruler].
This variable used in fixed format mode.")

(defvar fortran-column-ruler-tab
  "0       810        20        30        40        5\
0        60        70\n\
\[   ]|  { |    |    |    |    |    |    |    |    \
\|    |    |    |    |}\n"
  "String displayed above current line by \\[fortran-column-ruler].
This variable used in TAB format mode.")

(defvar fortran-mode-syntax-table nil
  "Syntax table in use in Fortran mode buffers.")

(defvar fortran-analyze-depth 100
  "Number of lines to scan to determine whether to use fixed or TAB \
format style.")

(defcustom fortran-break-before-delimiters t
  "*Non-nil causes filling to break lines before delimiters."
  :type 'boolean
  :group 'fortran)

(if fortran-mode-syntax-table
    ()
  (setq fortran-mode-syntax-table (make-syntax-table))
  ;; We might like `;' to be punctuation (g77 multi-statement lines),
  ;; but that screws abbrevs.
  (modify-syntax-entry ?\; "w" fortran-mode-syntax-table)
  (modify-syntax-entry ?\r " " fortran-mode-syntax-table)
  (modify-syntax-entry ?+ "." fortran-mode-syntax-table)
  (modify-syntax-entry ?- "." fortran-mode-syntax-table)
  (modify-syntax-entry ?= "." fortran-mode-syntax-table)
  (modify-syntax-entry ?* "." fortran-mode-syntax-table)
  (modify-syntax-entry ?/ "." fortran-mode-syntax-table)
  (modify-syntax-entry ?\' "\"" fortran-mode-syntax-table)
  (modify-syntax-entry ?\" "\"" fortran-mode-syntax-table)
  (modify-syntax-entry ?\\ "/" fortran-mode-syntax-table)
  ;; This might be better as punctuation, as for C, but this way you
  ;; can treat floating-point numbers as symbols.
  (modify-syntax-entry ?. "_" fortran-mode-syntax-table) ; e.g. `a /= b'
  (modify-syntax-entry ?_ "_" fortran-mode-syntax-table)
  (modify-syntax-entry ?$ "_" fortran-mode-syntax-table) ; esp. VMSisms
  (modify-syntax-entry ?\! "<" fortran-mode-syntax-table)
  (modify-syntax-entry ?\n ">" fortran-mode-syntax-table))

;; Comments are real pain in Fortran because there is no way to represent the
;; standard comment syntax in an Emacs syntax table (we can for VAX-style).
;; Therefore an unmatched quote in a standard comment will throw fontification
;; off on the wrong track.  So we do syntactic fontification with regexps.

;; Regexps done by simon@gnu with help from Ulrik Dickow <dickow@nbi.dk> and
;; probably others Si's forgotten about (sorry).

(defconst fortran-font-lock-keywords-1 nil
  "Subdued level highlighting for Fortran mode.")

(defconst fortran-font-lock-keywords-2 nil
  "Medium level highlighting for Fortran mode.")

(defconst fortran-font-lock-keywords-3 nil
  "Gaudy level highlighting for Fortran mode.")

(defun fortran-fontify-string (limit)
  (let ((match (match-string 1)))
    (cond ((string= "'" match)
	   (re-search-forward "\\([^'\n]*'?\\)" limit))
	  ((string= "\"" match)
	   (re-search-forward "\\([^\"\n]*\"?\\)" limit)))))

(let ((comment-chars "c!*d")		; `d' for `debugging' comments
      (fortran-type-types
       (eval-when-compile
	 (let ((re (regexp-opt
		    (let ((simple-types
			   '("character" "byte" "integer" "logical"
			     "none" "real" "complex"
			     "double precision" "double complex"))
			  (structured-types '("structure" "union" "map"))
			  (other-types '("record" "dimension"
					 "parameter" "common" "save"
					 "external" "intrinsic" "data"
					 "equivalence")))
		      (append
		       (mapcar (lambda (x) (concat "implicit " x))
			       simple-types)
		       simple-types
		       (mapcar (lambda (x) (concat "end " x))
			       structured-types)
		       structured-types
		       other-types)))))
	   ;; In the optimized regexp above, replace spaces by regexp
	   ;; for optional whitespace, which regexp-opt would have
	   ;; escaped.
	   (mapconcat #'identity (split-string re) "[ \t]*"))))
      (fortran-keywords
       (eval-when-compile
         (regexp-opt '("continue" "format" "end" "enddo" "if" "then"
                       "else" "endif" "elseif" "while" "inquire" "stop"
                       "return" "include" "open" "close" "read" "write"
                       "format" "print" "select" "case" "cycle" "exit"))))
      (fortran-logicals
       (eval-when-compile
         (regexp-opt '("and" "or" "not" "lt" "le" "eq" "ge" "gt" "ne"
                       "true" "false")))))

  (setq fortran-font-lock-keywords-1
        (list
         ;;
         ;; Fontify syntactically (assuming strings cannot be quoted
         ;; or span lines).
         (cons (concat "^[" comment-chars "].*") 'font-lock-comment-face)
         '(fortran-match-!-comment . font-lock-comment-face)
         (list (concat "^[^" comment-chars "\t\n]" (make-string 71 ?.)
                       "\\(.*\\)")
               '(1 font-lock-comment-face))
  	 '("\\(\\s\"\\)"		; single- or double-quoted string
  	   (1 font-lock-string-face)
  	   (fortran-fontify-string nil nil (1 font-lock-string-face)))
         ;;
         ;; Program, subroutine and function declarations, plus calls.
         (list (concat "\\<\\(block[ \t]*data\\|call\\|entry\\|function\\|"
                       "program\\|subroutine\\)\\>[ \t]*\\(\\sw+\\)?")
               '(1 font-lock-keyword-face)
               '(2 font-lock-function-name-face nil t))))

  (setq fortran-font-lock-keywords-2
        (append fortran-font-lock-keywords-1
                (list
                 ;;
                 ;; Fontify all type specifiers (must be first; see below).
                 (cons (concat "\\<\\(" fortran-type-types "\\)\\>")
                       'font-lock-type-face)
                 ;;
                 ;; Fontify all builtin keywords (except logical, do
                 ;; and goto; see below).
                 (concat "\\<\\(" fortran-keywords "\\)\\>")
                 ;;
                 ;; Fontify all builtin operators.
                 (concat "\\.\\(" fortran-logicals "\\)\\.")
                 ;;
                 ;; Fontify do/goto keywords and targets, and goto tags.
                 (list "\\<\\(do\\|go *to\\)\\>[ \t]*\\([0-9]+\\)?"
                       '(1 font-lock-keyword-face)
                       '(2 font-lock-constant-face nil t))
                 (cons "^ *\\([0-9]+\\)" 'font-lock-constant-face))))

  (setq fortran-font-lock-keywords-3
        (append
         ;;
         ;; The list `fortran-font-lock-keywords-1'.
         fortran-font-lock-keywords-1
         ;;
         ;; Fontify all type specifiers plus their declared items.
         (list
          (list (concat "\\<\\(" fortran-type-types "\\)\\>[ \t(/]*\\(*\\)?")
                ;; Fontify the type specifier.
                '(1 font-lock-type-face)
                ;; Fontify each declaration item (or just the /.../ block name).
                `(font-lock-match-c-style-declaration-item-and-skip-to-next
                  ;; Start after any *(...) expression.
                  (condition-case nil
		      (and (and (match-beginning ,(+ 2 (regexp-opt-depth
							fortran-type-types)))
				(forward-sexp))
			   (forward-sexp))
		    (error nil))
                  ;; No need to clean up.
                  nil
                  ;; Fontify as a variable name, functions are
                  ;; fontified elsewhere.
                  (1 font-lock-variable-name-face nil t))))
         ;;
         ;; Things extra to `fortran-font-lock-keywords-3'
         ;; (must be done first).
         (list
          ;;
          ;; Fontify goto-like `err=label'/`end=label' in read/write
          ;; statements.
          '(", *\\(e\\(nd\\|rr\\)\\)\\> *\\(= *\\([0-9]+\\)\\)?"
            (1 font-lock-keyword-face) (4 font-lock-constant-face nil t))
          ;;
          ;; Highlight standard continuation character and in a
          ;; TAB-formatted line.
          '("^     \\([^ 0]\\)" 1 font-lock-string-face)
          '("^\t\\([1-9]\\)" 1 font-lock-string-face))
         ;;
         ;; The list `fortran-font-lock-keywords-2' less that for types
         ;; (see above).
         (cdr (nthcdr (length fortran-font-lock-keywords-1)
                      fortran-font-lock-keywords-2)))))

(defvar fortran-font-lock-keywords fortran-font-lock-keywords-1
  "Default expressions to highlight in Fortran mode.")

(defvar fortran-imenu-generic-expression
  ;; These patterns could be confused by sequence nos. in cols 72+ and
  ;; don't allow continuations everywhere.
  (list
   (list
    nil
    ;; Lines below are: 1. leading whitespace; 2. function
    ;; declaration with optional type, e.g. `real', `real*4',
    ;; character(*), `double precision' and possible statement
    ;; continuation; 3. untyped declarations; 4. the variable to
    ;; index.  [This will be fooled by `end function' allowed by G77.
    ;; Also, it assumes sensible whitespace is employed.]
    (concat "^\\s-+\\(\
\\(\\sw\\|\\s-\\|[*()+]\\)*\
\\<function\\|subroutine\\|entry\\|block\\s-*data\\|program\\)\
[ \t" fortran-continuation-string "]+\
\\(\\sw+\\)")
    3)
   ;; Un-named block data
   (list nil "^\\s-+\\(block\\s-*data\\)\\s-*$" 1))
  "Imenu generic expression for `imenu-default-create-index-function'.")

(defvar fortran-mode-map ()
  "Keymap used in Fortran mode.")
(if fortran-mode-map
    ()
  (setq fortran-mode-map (make-sparse-keymap))
  (define-key fortran-mode-map ";" 'fortran-abbrev-start)
  (define-key fortran-mode-map "\C-c;" 'fortran-comment-region)
  (define-key fortran-mode-map "\M-\C-a" 'beginning-of-fortran-subprogram)
  (define-key fortran-mode-map "\M-\C-e" 'end-of-fortran-subprogram)
  (define-key fortran-mode-map "\M-;" 'fortran-indent-comment)
  (define-key fortran-mode-map "\M-\C-h" 'mark-fortran-subprogram)
  (define-key fortran-mode-map "\M-\n" 'fortran-split-line)
  (define-key fortran-mode-map "\n" 'fortran-indent-new-line)
  (define-key fortran-mode-map "\M-\C-q" 'fortran-indent-subprogram)
  (define-key fortran-mode-map "\C-c\C-w" 'fortran-window-create-momentarily)
  (define-key fortran-mode-map "\C-c\C-r" 'fortran-column-ruler)
  (define-key fortran-mode-map "\C-c\C-p" 'fortran-previous-statement)
  (define-key fortran-mode-map "\C-c\C-n" 'fortran-next-statement)
  (define-key fortran-mode-map "\C-c\C-d" 'fortran-join-line) ; like f90
  (define-key fortran-mode-map "\M-^" 'fortran-join-line) ; subvert delete-indentation
  (define-key fortran-mode-map "\C-xnd" 'fortran-narrow-to-subprogram)
  ;(define-key fortran-mode-map "\t" 'fortran-indent-line)
  (define-key fortran-mode-map "0" 'fortran-electric-line-number)
  (define-key fortran-mode-map "1" 'fortran-electric-line-number)
  (define-key fortran-mode-map "2" 'fortran-electric-line-number)
  (define-key fortran-mode-map "3" 'fortran-electric-line-number)
  (define-key fortran-mode-map "4" 'fortran-electric-line-number)
  (define-key fortran-mode-map "5" 'fortran-electric-line-number)
  (define-key fortran-mode-map "6" 'fortran-electric-line-number)
  (define-key fortran-mode-map "7" 'fortran-electric-line-number)
  (define-key fortran-mode-map "8" 'fortran-electric-line-number)
  (define-key fortran-mode-map "9" 'fortran-electric-line-number)

  ;; Menu
  (unless (boundp 'fortran-mode-menu)
    (easy-menu-define
     fortran-mode-menu fortran-mode-map ""
     '("Fortran"
       ["Toggle Auto-fill" fortran-auto-fill-mode :style toggle
        :selected (eq auto-fill-function 'fortran-do-auto-fill)]
       ["Toggle abbrev-mode" abbrev-mode :style toggle :selected abbrev-mode]
       "----"
       ["Comment-out Region" fortran-comment-region mark-active]
       ["Uncomment-out region"
        (fortran-comment-region (region-beginning) (region-end) 1)
        mark-active]
       ["Indent Region" indent-region mark-active]
       ["Indent Subprogram" fortran-indent-subprogram t]
       "----"
       ["Beginning of Subprogram" beginning-of-fortran-subprogram t]
       ["End of Subprogram" end-of-fortran-subprogram t]
       ("Mark"
        ["Subprogram" mark-fortran-subprogram t]
        ["IF Block" fortran-mark-if t]
        ["DO Block" fortran-mark-do t])
       ["Narrow to Subprogram" fortran-narrow-to-subprogram t]
       ["Widen" widen t]
       "----"
       ["Temporary column ruler" fortran-column-ruler t]
       ["72-column window" fortran-window-create t]
       ["Full Width Window"
        (enlarge-window-horizontally (- (frame-width) (window-width)))
        (< (window-width) (frame-width))]
       ["Momentary 72-column window" fortran-window-create-momentarily t]
       "----"
       ["Break Line at Point" fortran-split-line t]
       ["Join Line" fortran-join-line t]
       ["Fill Statement/Comment" fill-paragraph  t]
       "----"
       ["Add imenu menu"
        (progn (imenu-add-menubar-index)
               ;; Prod menu bar to update -- is this the right way?
               (menu-bar-mode 1))
        (not (and (boundp 'imenu--index-alist)
		  imenu--index-alist))]))))

(defvar fortran-mode-abbrev-table nil)
(if fortran-mode-abbrev-table
    ()
  (let ((ac abbrevs-changed))
    (define-abbrev-table 'fortran-mode-abbrev-table ())
    (define-abbrev fortran-mode-abbrev-table  ";au"  "automatic" nil)
    (define-abbrev fortran-mode-abbrev-table  ";b"   "byte" nil)
    (define-abbrev fortran-mode-abbrev-table  ";bd"  "block data" nil)
    (define-abbrev fortran-mode-abbrev-table  ";ch"  "character" nil)
    (define-abbrev fortran-mode-abbrev-table  ";cl"  "close" nil)
    (define-abbrev fortran-mode-abbrev-table  ";c"   "continue" nil)
    (define-abbrev fortran-mode-abbrev-table  ";cm"  "common" nil)
    (define-abbrev fortran-mode-abbrev-table  ";cx"  "complex" nil)
    (define-abbrev fortran-mode-abbrev-table  ";df"  "define" nil)
    (define-abbrev fortran-mode-abbrev-table  ";di"  "dimension" nil)
    (define-abbrev fortran-mode-abbrev-table  ";do"  "double" nil)
    (define-abbrev fortran-mode-abbrev-table  ";dc"  "double complex" nil)
    (define-abbrev fortran-mode-abbrev-table  ";dp"  "double precision" nil)
    (define-abbrev fortran-mode-abbrev-table  ";dw"  "do while" nil)
    (define-abbrev fortran-mode-abbrev-table  ";e"   "else" nil)
    (define-abbrev fortran-mode-abbrev-table  ";ed"  "enddo" nil)
    (define-abbrev fortran-mode-abbrev-table  ";el"  "elseif" nil)
    (define-abbrev fortran-mode-abbrev-table  ";en"  "endif" nil)
    (define-abbrev fortran-mode-abbrev-table  ";eq"  "equivalence" nil)
    (define-abbrev fortran-mode-abbrev-table  ";ew"  "endwhere" nil)
    (define-abbrev fortran-mode-abbrev-table  ";ex"  "external" nil)
    (define-abbrev fortran-mode-abbrev-table  ";ey"  "entry" nil)
    (define-abbrev fortran-mode-abbrev-table  ";f"   "format" nil)
    (define-abbrev fortran-mode-abbrev-table  ";fa"  ".false." nil)
    (define-abbrev fortran-mode-abbrev-table  ";fu"  "function" nil)
    (define-abbrev fortran-mode-abbrev-table  ";g"   "goto" nil)
    (define-abbrev fortran-mode-abbrev-table  ";im"  "implicit" nil)
    (define-abbrev fortran-mode-abbrev-table  ";ib"  "implicit byte" nil)
    (define-abbrev fortran-mode-abbrev-table  ";ic"  "implicit complex" nil)
    (define-abbrev fortran-mode-abbrev-table  ";ich" "implicit character" nil)
    (define-abbrev fortran-mode-abbrev-table  ";ii"  "implicit integer" nil)
    (define-abbrev fortran-mode-abbrev-table  ";il"  "implicit logical" nil)
    (define-abbrev fortran-mode-abbrev-table  ";ir"  "implicit real" nil)
    (define-abbrev fortran-mode-abbrev-table  ";inc" "include" nil)
    (define-abbrev fortran-mode-abbrev-table  ";in"  "integer" nil)
    (define-abbrev fortran-mode-abbrev-table  ";intr" "intrinsic" nil)
    (define-abbrev fortran-mode-abbrev-table  ";l"   "logical" nil)
    (define-abbrev fortran-mode-abbrev-table  ";n"   "namelist" nil)
    (define-abbrev fortran-mode-abbrev-table  ";o"   "open" nil) ; was ;op
    (define-abbrev fortran-mode-abbrev-table  ";pa"  "parameter" nil)
    (define-abbrev fortran-mode-abbrev-table  ";pr"  "program" nil)
    (define-abbrev fortran-mode-abbrev-table  ";ps"  "pause" nil)
    (define-abbrev fortran-mode-abbrev-table  ";p"   "print" nil)
    (define-abbrev fortran-mode-abbrev-table  ";rc"  "record" nil)
    (define-abbrev fortran-mode-abbrev-table  ";re"  "real" nil)
    (define-abbrev fortran-mode-abbrev-table  ";r"   "read" nil)
    (define-abbrev fortran-mode-abbrev-table  ";rt"  "return" nil)
    (define-abbrev fortran-mode-abbrev-table  ";rw"  "rewind" nil)
    (define-abbrev fortran-mode-abbrev-table  ";s"   "stop" nil)
    (define-abbrev fortran-mode-abbrev-table  ";sa"  "save" nil)
    (define-abbrev fortran-mode-abbrev-table  ";st"  "structure" nil)
    (define-abbrev fortran-mode-abbrev-table  ";sc"  "static" nil)
    (define-abbrev fortran-mode-abbrev-table  ";su"  "subroutine" nil)
    (define-abbrev fortran-mode-abbrev-table  ";tr"  ".true." nil)
    (define-abbrev fortran-mode-abbrev-table  ";ty"  "type" nil)
    (define-abbrev fortran-mode-abbrev-table  ";vo"  "volatile" nil)
    (define-abbrev fortran-mode-abbrev-table  ";w"   "write" nil)
    (define-abbrev fortran-mode-abbrev-table  ";wh"  "where" nil)
    (setq abbrevs-changed ac)))

(eval-when-compile			; silence compiler
  (defvar imenu-case-fold-search)
  (defvar imenu-syntax-alist))

;;;###autoload
(defun fortran-mode ()
  "Major mode for editing Fortran code.
\\[fortran-indent-line] indents the current Fortran line correctly.
DO statements must not share a common CONTINUE.

Type ;? or ;\\[help-command] to display a list of built-in abbrevs for
Fortran keywords.

Key definitions:
\\{fortran-mode-map}

Variables controlling indentation style and extra features:

 `comment-start'
    Normally nil in Fortran mode.  If you want to use comments
    starting with `!', set this to the string \"!\".
 `fortran-do-indent'
    Extra indentation within do blocks.  (default 3)
 `fortran-if-indent'
    Extra indentation within if blocks.  (default 3)
 `fortran-structure-indent'
    Extra indentation within structure, union, map and interface blocks.
    (default 3)
 `fortran-continuation-indent'
    Extra indentation applied to continuation statements.  (default 5)
 `fortran-comment-line-extra-indent'
    Amount of extra indentation for text within full-line comments.  (default 0)
 `fortran-comment-indent-style'
    nil    means don't change indentation of text in full-line comments,
    fixed  means indent that text at `fortran-comment-line-extra-indent' beyond
           the value of `fortran-minimum-statement-indent-fixed' (for fixed
           format continuation style) or `fortran-minimum-statement-indent-tab'
           (for TAB format continuation style).
    relative  means indent at `fortran-comment-line-extra-indent' beyond the
 	      indentation for a line of code.
    (default 'fixed)
 `fortran-comment-indent-char'
    Single-character string to be inserted instead of space for
    full-line comment indentation.  (default \" \")
 `fortran-minimum-statement-indent-fixed'
    Minimum indentation for Fortran statements in fixed format mode.  (def.6)
 `fortran-minimum-statement-indent-tab'
    Minimum indentation for Fortran statements in TAB format mode.  (default 9)
 `fortran-line-number-indent'
    Maximum indentation for line numbers.  A line number will get
    less than this much indentation if necessary to avoid reaching
    column 5.  (default 1)
 `fortran-check-all-num-for-matching-do'
    Non-nil causes all numbered lines to be treated as possible \"continue\"
    statements.  (default nil)
 `fortran-blink-matching-if'
    Non-nil causes \\[fortran-indent-line] on an ENDIF statement to blink on
    matching IF.  Also, from an ENDDO statement, blink on matching DO [WHILE]
    statement.  (default nil)
 `fortran-continuation-string'
    Single-character string to be inserted in column 5 of a continuation
    line.  (default \"$\")
 `fortran-comment-region'
    String inserted by \\[fortran-comment-region] at start of each line in
    region.  (default \"c$$$\")
 `fortran-electric-line-number'
    Non-nil causes line number digits to be moved to the correct column
    as typed.  (default t)
 `fortran-break-before-delimiters'
    Non-nil causes `fortran-fill' to break lines before delimiters.
    (default t)

Turning on Fortran mode calls the value of the variable `fortran-mode-hook'
with no args, if that value is non-nil."
  (interactive)
  (kill-all-local-variables)
  (setq local-abbrev-table fortran-mode-abbrev-table)
  (set-syntax-table fortran-mode-syntax-table)
  ;; Font Lock mode support.
  (make-local-variable 'font-lock-defaults)
  (setq font-lock-defaults '((fortran-font-lock-keywords
			      fortran-font-lock-keywords-1
			      fortran-font-lock-keywords-2
			      fortran-font-lock-keywords-3)
			     t t ((?/ . "$/") ("_$" . "w"))))
  (make-local-variable 'fortran-break-before-delimiters)
  (setq fortran-break-before-delimiters t)
  (make-local-variable 'indent-line-function)
  (setq indent-line-function 'fortran-indent-line)
  (make-local-variable 'comment-indent-function)
  (setq comment-indent-function 'fortran-comment-hook)
  (make-local-variable 'comment-line-start-skip)
  (setq comment-line-start-skip
	"^[Cc*]\\(\\([^ \t\n]\\)\\2\\2*\\)?[ \t]*\\|^#.*")
  (make-local-variable 'comment-line-start)
  (setq comment-line-start "c")
  (make-local-variable 'comment-start-skip)
  (setq comment-start-skip "![ \t]*")
  (make-local-variable 'comment-start)
  (setq comment-start nil)
  (make-local-variable 'require-final-newline)
  (setq require-final-newline t)
  (make-local-variable 'abbrev-all-caps)
  (setq abbrev-all-caps t)
  (make-local-variable 'indent-tabs-mode)
  (setq indent-tabs-mode nil)
;;;(setq abbrev-mode t) ; ?? (abbrev-mode 1) instead??
  (set (make-local-variable 'fill-column) 71)
  (use-local-map fortran-mode-map)
  (setq mode-name "Fortran")
  (setq major-mode 'fortran-mode)
  (make-local-variable 'fortran-comment-line-extra-indent)
  (make-local-variable 'fortran-minimum-statement-indent-fixed)
  (make-local-variable 'fortran-minimum-statement-indent-tab)
  (make-local-variable 'fortran-column-ruler-fixed)
  (make-local-variable 'fortran-column-ruler-tab)
  (setq fortran-tab-mode-string " TAB-format")
  (setq indent-tabs-mode (fortran-analyze-file-format))
  (setq imenu-case-fold-search t)
  (make-local-variable 'imenu-generic-expression)
  (setq imenu-generic-expression fortran-imenu-generic-expression)
  (setq imenu-syntax-alist '(("_$" . "w")))
  (set (make-local-variable 'fill-paragraph-function) 'fortran-fill-paragraph)
  (set (make-local-variable 'indent-line-function) 'fortran-indent-line)
  (set (make-local-variable 'indent-region-function)
       (lambda (start end)
         (let (fortran-blink-matching-if ; avoid blinking delay
               indent-region-function)
           (indent-region start end nil))))
  (run-hooks 'fortran-mode-hook))

(defun fortran-comment-hook ()
  (save-excursion
    (skip-chars-backward " \t")
    (max (+ 1 (current-column))
	 comment-column)))

(defun fortran-indent-comment ()
  "Align or create comment on current line.
Existing comments of all types are recognized and aligned.
If the line has no comment, a side-by-side comment is inserted and aligned
if the value of  `comment-start'  is not nil.
Otherwise, a separate-line comment is inserted, on this line
or on a new line inserted before this line if this line is not blank."
  (interactive)
  (beginning-of-line)
  ;; Recognize existing comments of either kind.
  (cond ((looking-at comment-line-start-skip)
	 (fortran-indent-line))
	((fortran-find-comment-start-skip) ; catches any inline comment and
					; leaves point after comment-start-skip
	 (if comment-start-skip
	     (progn (goto-char (match-beginning 0))
		    (if (not (= (current-column) (fortran-comment-hook)))
			(progn (delete-horizontal-space)
			       (indent-to (fortran-comment-hook)))))
	   (end-of-line)))        ; otherwise goto end of line or sth else?
	;; No existing comment.
	;; If side-by-side comments are defined, insert one,
	;; unless line is now blank.
	((and comment-start (not (looking-at "^[ \t]*$")))
	 (end-of-line)
	 (delete-horizontal-space)
	 (indent-to (fortran-comment-hook))
	 (insert comment-start))
	;; Else insert separate-line comment, making a new line if nec.
	(t
	 (if (looking-at "^[ \t]*$")
	     (delete-horizontal-space)
	   (beginning-of-line)
	   (insert "\n")
	   (forward-char -1))
	 (insert comment-line-start)
	 (insert-char (if (stringp fortran-comment-indent-char)
			  (aref fortran-comment-indent-char 0)
			fortran-comment-indent-char)
		      (- (fortran-calculate-indent) (current-column))))))

(defun fortran-comment-region (beg-region end-region arg)
  "Comments every line in the region.
Puts `fortran-comment-region' at the beginning of every line in the region.
BEG-REGION and END-REGION are args which specify the region boundaries.
With non-nil ARG, uncomments the region."
  (interactive "*r\nP")
  (let ((end-region-mark (make-marker)) (save-point (point-marker)))
    (set-marker end-region-mark end-region)
    (goto-char beg-region)
    (beginning-of-line)
    (if (not arg)			;comment the region
	(progn (insert fortran-comment-region)
	       (while (and  (= (forward-line 1) 0)
			    (< (point) end-region-mark))
		 (insert fortran-comment-region)))
      (let ((com (regexp-quote fortran-comment-region))) ;uncomment the region
	(if (looking-at com)
	    (delete-region (point) (match-end 0)))
	(while (and  (= (forward-line 1) 0)
		     (< (point) end-region-mark))
	  (if (looking-at com)
	      (delete-region (point) (match-end 0))))))
    (goto-char save-point)
    (set-marker end-region-mark nil)
    (set-marker save-point nil)))

(defun fortran-abbrev-start ()
  "Typing ;\\[help-command] or ;? lists all the Fortran abbrevs.
Any other key combination is executed normally."
  (interactive)
  (let (c)
    (insert last-command-char)
    (if (or (eq (setq c (read-event)) ??)    ;insert char if not equal to `?'
	    (eq c help-char))
	(fortran-abbrev-help)
      (setq unread-command-events (list c)))))

(defun fortran-abbrev-help ()
  "List the currently defined abbrevs in Fortran mode."
  (interactive)
  (message "Listing abbrev table...")
  (display-buffer (fortran-prepare-abbrev-list-buffer))
  (message "Listing abbrev table...done"))

(defun fortran-prepare-abbrev-list-buffer ()
  (save-excursion
    (set-buffer (get-buffer-create "*Abbrevs*"))
    (erase-buffer)
    (insert-abbrev-table-description 'fortran-mode-abbrev-table t)
    (goto-char (point-min))
    (set-buffer-modified-p nil)
    (edit-abbrevs-mode))
  (get-buffer-create "*Abbrevs*"))

(defun fortran-column-ruler ()
  "Insert a column ruler momentarily above current line, till next keystroke.
The ruler is defined by the value of `fortran-column-ruler-fixed' when in fixed
format mode, and `fortran-column-ruler-tab' when in TAB format mode.
The key typed is executed unless it is SPC."
  (interactive)
  (momentary-string-display
   (if indent-tabs-mode
       fortran-column-ruler-tab
     fortran-column-ruler-fixed)
   (save-excursion
     (beginning-of-line)
     (if (eq (window-start (selected-window))
	     (window-point (selected-window)))
	 (progn (forward-line) (point))
       (point)))
   nil "Type SPC or any command to erase ruler."))

(defun fortran-window-create ()
  "Make the window 72 columns wide.
See also `fortran-window-create-momentarily'."
  (interactive)
  (condition-case error
      (progn
	(let ((window-min-width 2))
	  (if (< (window-width) (frame-width))
	      (enlarge-window-horizontally (- (frame-width)
					      (window-width) 1)))
	  (let* ((window-edges (window-edges))
		 (scroll-bar-width (- (nth 2 window-edges)
				      (car window-edges)
				      (window-width))))
	    (split-window-horizontally (+ 72 scroll-bar-width)))
	  (other-window 1)
	  (switch-to-buffer " fortran-window-extra" t)
	  (select-window (previous-window))))
    (error (message "No room for Fortran window.")
	   'error)))

(defun fortran-window-create-momentarily (&optional arg)
  "Momentarily make the window 72 columns wide.
Optional ARG non-nil and non-unity disables the momentary feature.
See also `fortran-window-create'."
  (interactive "p")
  (if (or (not arg)
	  (= arg 1))
      (save-window-excursion
	(if (not (equal (fortran-window-create) 'error))
	    (progn (message "Type SPC to continue editing.")
		   (let ((char (read-event)))
		     (or (equal char (string-to-char " "))
			 (setq unread-command-events (list char)))))))
    (fortran-window-create)))

(defun fortran-split-line ()
  "Break line at point and insert continuation marker and alignment."
  (interactive)
  (delete-horizontal-space)
  (if (save-excursion (beginning-of-line) (looking-at comment-line-start-skip))
      (insert "\n" comment-line-start " ")
    (if indent-tabs-mode
	(insert "\n\t" (fortran-numerical-continuation-char))
      (insert "\n " fortran-continuation-string))) ; Space after \n important
  (fortran-indent-line))		; when the cont string is C, c or *.

(defun fortran-remove-continuation ()
  (if (looking-at "\\(     [^ 0\n]\\|\t[1-9]\\|&\\)")
      (progn (replace-match "")
	     (delete-indentation)
	     t)))

(defun fortran-join-line (arg)
  "Join current line to the previous one and re-indent.
With a prefix argument, repeat this operation that many times.
If the prefix argument ARG is negative, join the next -ARG lines.
Continuation lines are correctly handled."
  (interactive "*p")
  (save-excursion
    (when (> 0 arg)
      (setq arg (- arg))
      (forward-line arg))
    (while (not (zerop arg))
      (beginning-of-line)
      (or (fortran-remove-continuation)
          (delete-indentation))
      (setq arg (1- arg)))
    (fortran-indent-line)))

(defun fortran-numerical-continuation-char ()
  "Return a digit for tab-digit style of continuation lines.
If, previous line is a tab-digit continuation line, returns that digit
plus one.  Otherwise return 1.  Zero not allowed."
  (save-excursion
    (forward-line -1)
    (if (looking-at "\t[1-9]")
	(+ ?1 (% (- (char-after (+ (point) 1)) ?0) 9))
      ?1)))

(defun delete-horizontal-regexp (chars)
  "Delete all characters in CHARS around point.
CHARS is like the inside of a [...] in a regular expression
except that ] is never special and \ quotes ^, - or \."
  (interactive "*s")
  (skip-chars-backward chars)
  (delete-region (point) (progn (skip-chars-forward chars) (point))))

(put 'fortran-electric-line-number 'delete-selection t)
(defun fortran-electric-line-number (arg)
  "Self insert, but if part of a Fortran line number indent it automatically.
Auto-indent does not happen if a numeric ARG is used."
  (interactive "P")
  (if (or arg (not fortran-electric-line-number))
      (if arg
	  (self-insert-command (prefix-numeric-value arg))
	(self-insert-command 1))
    (if (or (and (= 5 (current-column))
		 (save-excursion
		   (beginning-of-line)
		   (looking-at "     ")));In col 5 with only spaces to left.
	    (and (= (if indent-tabs-mode
			fortran-minimum-statement-indent-tab
		      fortran-minimum-statement-indent-fixed) (current-column))
		 (save-excursion
		   (beginning-of-line)
		   (looking-at "\t"));In col 8 with a single tab to the left.
		 (not (or (eq last-command 'fortran-indent-line)
			  (eq last-command
			      'fortran-indent-new-line))))
	    (save-excursion
	      (re-search-backward "[^ \t0-9]"
				  (save-excursion
				    (beginning-of-line)
				    (point))
				  t))	;not a line number
	    (looking-at "[0-9]"))	;within a line number
	(self-insert-command (prefix-numeric-value arg))
      (skip-chars-backward " \t")
      (insert last-command-char)
      (fortran-indent-line))))

(defvar fortran-end-prog-re1
  "\\(end\\|contains\\)\
\\([ \t]*\\(program\\|subroutine\\|function\\|block[ \t]*data\\)\\>\
\\([ \t]*\\(\\sw\\|\\s_\\)+\\)?\\)?"
"Regexp possibly marking subprogram end.")
(defvar fortran-end-prog-re
  (concat "^[ \t0-9]*" fortran-end-prog-re1)
  "Regexp possibly marking subprogram end.")

(defun fortran-check-end-prog-re ()
  "Check a preliminary match against `fortran-end-prog-re'."
  ;; Having got a possible match for the subprogram end, we need a
  ;; match of whitespace, avoiding possible column 73+ stuff.
  (save-match-data
    (string-match "^\\s-*\\(\\'\\|\\s<\\)"
		  (buffer-substring (match-end 0)
				    (min (line-end-position)
					 (+ 72 (line-beginning-position)))))))

;; Note that you can't just check backwards for `subroutine' &c in
;; case of un-marked main programs not at the start of the file.
;; Nov 2013 MvS updated to jump past "end interface"
(defun beginning-of-fortran-subprogram ()
  "Moves point to the beginning of the current Fortran subprogram."
  (interactive)
  (let ((case-fold-search t))
    (beginning-of-line -1)
    (if (catch 'ok
	  (while (re-search-backward fortran-end-prog-re nil 'move)
	    (if (fortran-check-end-prog-re)
		(throw 'ok t))))
	(forward-line))
    (if (looking-at   "[ \t]+end[ \t]+interface") 
      (beginning-of-fortran-subprogram))))

;; Nov 2013 MvS updated to jump past "end interface"
(defun end-of-fortran-subprogram ()
  "Moves point to the end of the current Fortran subprogram."
  (interactive)
  (let ((case-fold-search t))
    (if (save-excursion			; already sitting on END
	  (beginning-of-line)
	  (and (looking-at fortran-end-prog-re)
	       (fortran-check-end-prog-re)))
	(forward-line)
      (beginning-of-line 2) ; move past END
      (catch 'ok
	(while (re-search-forward fortran-end-prog-re nil 'move)
	  (if (fortran-check-end-prog-re)
	      (throw 'ok t))))
      (goto-char (match-beginning 0))
      (forward-line))
    (if (looking-at   "[ \t]+end[ \t]+interface") 
      (end-of-fortran-subprogram)))
  (point)
)

(defun mark-fortran-subprogram ()
  "Put mark at end of Fortran subprogram, point at beginning.
The marks are pushed."
  (interactive)
  (end-of-fortran-subprogram)
  (push-mark (point) nil t)
  (beginning-of-fortran-subprogram))

(defun fortran-previous-statement ()
  "Moves point to beginning of the previous Fortran statement.
Returns `first-statement' if that statement is the first
non-comment Fortran statement in the file, and nil otherwise."
  (interactive)
  (let (not-first-statement continue-test)
    (beginning-of-line)
    (setq continue-test
	  (and
	   (not (looking-at comment-line-start-skip))
	   (or (looking-at
	        (concat "[ \t]*" (regexp-quote fortran-continuation-string)))
	       (or (looking-at "     [^ 0\n]")
		   (looking-at "\t[1-9]")))))
    (while (and (setq not-first-statement (= (forward-line -1) 0))
		(or (looking-at comment-line-start-skip)
		    (looking-at "[ \t]*$")
		    (looking-at "     [^ 0\n]")
		    (looking-at "\t[1-9]")
		    (looking-at (concat "[ \t]*"  comment-start-skip)))))
    (cond ((and continue-test
		(not not-first-statement))
	   (message "Incomplete continuation statement."))
	  (continue-test
	   (fortran-previous-statement))
	  ((not not-first-statement)
	   'first-statement))))

(defun fortran-next-statement ()
  "Moves point to beginning of the next Fortran statement.
Returns `last-statement' if that statement is the last
non-comment Fortran statement in the file, and nil otherwise."
  (interactive)
  (let (not-last-statement)
    (beginning-of-line)
    (while (and (setq not-last-statement
		      (and (= (forward-line 1) 0)
			   (not (eobp))))
 		(or (looking-at comment-line-start-skip)
 		    (looking-at "[ \t]*$")
 		    (looking-at "     [^ 0\n]")
 		    (looking-at "\t[1-9]")
 		    (looking-at (concat "[ \t]*"  comment-start-skip)))))
    (if (not not-last-statement)
 	'last-statement)))

(defun fortran-narrow-to-subprogram ()
  "Make text outside the current subprogram invisible.
The subprogram visible is the one that contains or follows point."
  (interactive)
  (save-excursion
    (mark-fortran-subprogram)
    (narrow-to-region (point) (mark))))

(defmacro fortran-with-subprogram-narrowing (&rest forms)
  "Execute FORMS with buffer temporarily narrowed to current subprogram.
Doesn't push a mark."
  `(save-restriction
     (save-excursion
       (narrow-to-region (progn
			   (beginning-of-fortran-subprogram)
			   (point))
			 (progn
			   (end-of-fortran-subprogram)
			   (point))))
     ,@forms))

(defun fortran-blink-matching-if ()
  "From an ENDIF statement, blink the matching IF statement."
  (let ((top-of-window (window-start))
	(endif-point (point))
	(case-fold-search t)
	matching-if
	message)
    (if (save-excursion (beginning-of-line)
			(skip-chars-forward " \t0-9")
			(looking-at "e\\(nd[ \t]*if\\|lse\\([ \t]*if\\)?\\)\\b"))
	(progn
          (if (not (setq matching-if (fortran-beginning-if)))
              (setq message "No matching if.")
            (if (< matching-if top-of-window)
                (save-excursion
                  (goto-char matching-if)
                  (beginning-of-line)
                  (setq message
                        (concat "Matches "
                                (buffer-substring
                                 (point) (progn (end-of-line) (point))))))))
	  (if message
	      (message "%s" message)
	    (goto-char matching-if)
	    (sit-for 1)
	    (goto-char endif-point))))))

(defun fortran-blink-matching-do ()
  "From an ENDDO statement, blink the matching DO or DO WHILE statement."
  ;; This is basically copied from fortran-blink-matching-if.
  (let ((top-of-window (window-start))
	(enddo-point (point))
	(case-fold-search t)
	matching-do
	message)
    (if (save-excursion (beginning-of-line)
			(skip-chars-forward " \t0-9")
			(looking-at "end[ \t]*do\\b"))
	(progn
          (if (not (setq matching-do (fortran-beginning-do)))
              (setq message "No matching do.")
            (if (< matching-do top-of-window)
                (save-excursion
                  (goto-char matching-do)
                  (beginning-of-line)
                  (setq message
                        (concat "Matches "
                                (buffer-substring
                                 (point) (progn (end-of-line) (point))))))))
	  (if message
	      (message "%s" message)
	    (goto-char matching-do)
	    (sit-for 1)
	    (goto-char enddo-point))))))

(defun fortran-mark-do ()
  "Put mark at end of Fortran DO [WHILE]-ENDDO construct, point at beginning.
The marks are pushed."
  (interactive)
  (let (enddo-point do-point)
    (if (setq enddo-point (fortran-end-do))
        (if (not (setq do-point (fortran-beginning-do)))
            (message "No matching do.")
          ;; Set mark, move point.
          (goto-char enddo-point)
          (push-mark)
          (goto-char do-point)))))

(defun fortran-end-do ()
  "Search forward for first unmatched ENDDO.
Return point or nil."
  (let ((case-fold-search t))
    (if (save-excursion (beginning-of-line)
			(skip-chars-forward " \t0-9")
			(looking-at "end[ \t]*do\\b"))
	;; Sitting on one.
	(match-beginning 0)
      ;; Search for one.
      (save-excursion
	(let ((count 1))
        (while (and (not (= count 0))
		      (not (eq (fortran-next-statement) 'last-statement))
		      ;; Keep local to subprogram
		      (not (and (looking-at fortran-end-prog-re)
				(fortran-check-end-prog-re))))

	    (skip-chars-forward " \t0-9")
	    (cond ((looking-at "end[ \t]*do\\b")
		   (setq count (1- count)))
		  ((looking-at "\\(\\(\\sw\\|\\s_\\)+:[ \t]*\\)?do[ \t]+[^0-9]")
                 (setq count (+ count 1)))))
        (and (= count 0)
	       ;; All pairs accounted for.
	       (point)))))))

(defun fortran-beginning-do ()
  "Search backwards for first unmatched DO [WHILE].
Return point or nil."
  (let ((case-fold-search t))
    (if (save-excursion (beginning-of-line)
			(skip-chars-forward " \t0-9")
			(looking-at "\\(\\(\\sw\\|\\s_\\)+:[ \t]*\\)?do[ \t]+"))
	;; Sitting on one.
	(match-beginning 0)
      ;; Search for one.
      (save-excursion
	(let ((count 1))
        (while (and (not (= count 0))
		      (not (eq (fortran-previous-statement) 'first-statement))
		      ;; Keep local to subprogram
		      (not (and (looking-at fortran-end-prog-re)
				(fortran-check-end-prog-re))))

	    (skip-chars-forward " \t0-9")
	    (cond ((looking-at "\\(\\(\\sw\\|\\s_\\)+:[ \t]*\\)?do[ \t]+[^0-9]")
		   (setq count (1- count)))
		  ((looking-at "end[ \t]*do\\b")
		   (setq count (1+ count)))))

        (and (= count 0)
	       ;; All pairs accounted for.
	       (point)))))))

(defun fortran-mark-if ()
  "Put mark at end of Fortran IF-ENDIF construct, point at beginning.
The marks are pushed."
  (interactive)
  (let (endif-point if-point)
    (if (setq endif-point (fortran-end-if))
        (if (not (setq if-point (fortran-beginning-if)))
            (message "No matching if.")
          ;; Set mark, move point.
          (goto-char endif-point)
          (push-mark)
          (goto-char if-point)))))

(defvar fortran-if-start-re "\\(\\(\\sw\\|\\s_\\)+:[ \t]*\\)?if[ \t]*(")

(defun fortran-end-if ()
  "Search forwards for first unmatched ENDIF.
Return point or nil."
  (let ((case-fold-search t))
    (if (save-excursion (beginning-of-line)
			(skip-chars-forward " \t0-9")
			(looking-at "end[ \t]*if\\b"))
	;; Sitting on one.
	(match-beginning 0)
      ;; Search for one.  The point has been already been moved to first
      ;; letter on line but this should not cause troubles.
      (save-excursion
	(let ((count 1))
        (while (and (not (= count 0))
		      (not (eq (fortran-next-statement) 'last-statement))
		      ;; Keep local to subprogram.
		      (not (and (looking-at fortran-end-prog-re)
				(fortran-check-end-prog-re))))

	    (skip-chars-forward " \t0-9")
	    (cond ((looking-at "end[ \t]*if\\b")
                 (setq count (- count 1)))

		  ((looking-at fortran-if-start-re)
		   (save-excursion
		     (if (or
			  (looking-at ".*)[ \t]*then\\b[ \t]*[^ \t(=a-z0-9]")
			  (let (then-test) ; Multi-line if-then.
			    (while
                              (and (= (forward-line 1) 0)
				     ;; Search forward for then.
				     (or (looking-at "     [^ 0\n]")
					 (looking-at "\t[1-9]"))
				     (not
				      (setq then-test
					    (looking-at
					     ".*then\\b[ \t]*[^ \t(=a-z0-9]")))))
			    then-test))
                       (setq count (+ count 1)))))))

        (and (= count 0)
	       ;; All pairs accounted for.
	       (point)))))))

(defun fortran-beginning-if ()
  "Search backwards for first unmatched IF-THEN.
Return point or nil."
  (let ((case-fold-search t))
    (if (save-excursion
	  ;; May be sitting on multi-line if-then statement, first move to
	  ;; beginning of current statement.  Note: `fortran-previous-statement'
	  ;; moves to previous statement *unless* current statement is first
	  ;; one.  Only move forward if not first-statement.
	  (if (not (eq (fortran-previous-statement) 'first-statement))
	      (fortran-next-statement))
	  (skip-chars-forward " \t0-9")
	  (and
	   (looking-at fortran-if-start-re)
	   (save-match-data
	     (or (looking-at ".*)[ \t]*then\\b[ \t]*[^ \t(=a-z0-9]")
		 ;; Multi-line if-then.
		 (let (then-test)
		   (while
                     (and (= (forward-line 1) 0)
			    ;; Search forward for then.
			    (or (looking-at "     [^ 0\n]")
				(looking-at "\t[1-9]"))
			    (not
			     (setq then-test
				   (looking-at
				    ".*then\\b[ \t]*[^ \t(=a-z0-9]")))))
		   then-test)))))
	;; Sitting on one.
	(match-beginning 0)
      ;; Search for one.
      (save-excursion
	(let ((count 1))
        (while (and (not (= count 0))
		      (not (eq (fortran-previous-statement) 'first-statement))
		      ;; Keep local to subprogram.
		      (not (and (looking-at fortran-end-prog-re)
				(fortran-check-end-prog-re))))

	    (skip-chars-forward " \t0-9")
	    (cond ((looking-at fortran-if-start-re)
		   (save-excursion
		     (if (or
			  (looking-at ".*)[ \t]*then\\b[ \t]*[^ \t(=a-z0-9]")
			  (let (then-test) ; Multi-line if-then.
			    (while
                              (and (= (forward-line 1) 0)
				     ;; Search forward for then.
				     (or (looking-at "     [^ 0\n]")
					 (looking-at "\t[1-9]"))
				     (not
				      (setq then-test
					    (looking-at
					     ".*then\\b[ \t]*[^ \t(=a-z0-9]")))))
			    then-test))
                       (setq count (- count 1)))))
		  ((looking-at "end[ \t]*if\\b")
                 (setq count (+ count 1)))))

        (and (= count 0)
	       ;; All pairs accounted for.
	       (point)))))))

(defun fortran-indent-line ()
  "Indent current Fortran line based on its contents and on previous lines."
  (interactive)
  (let ((cfi (fortran-calculate-indent)))
    (save-excursion
      (beginning-of-line)
      (if (or (not (= cfi (fortran-current-line-indentation)))
	      (and (re-search-forward "^[ \t]*[0-9]+" (+ (point) 4) t)
		   (not (fortran-line-number-indented-correctly-p))))
	  (fortran-indent-to-column cfi)
	(beginning-of-line)
	(if (and (not (looking-at comment-line-start-skip))
		 (fortran-find-comment-start-skip))
	    (fortran-indent-comment))))
    ;; Never leave point in left margin.
    (if (< (current-column) cfi)
	(move-to-column cfi))
    (if (and auto-fill-function
	     (> (save-excursion (end-of-line) (current-column)) fill-column))
	(save-excursion
	  (end-of-line)
	  (fortran-fill)))
    (if fortran-blink-matching-if
        (progn
	  (fortran-blink-matching-if)
	  (fortran-blink-matching-do)))))

(defun fortran-indent-new-line ()
  "Reindent the current Fortran line, insert a newline and indent the newline.
An abbrev before point is expanded if variable `abbrev-mode' is non-nil."
  (interactive)
  (if abbrev-mode (expand-abbrev))
  (save-excursion
    (beginning-of-line)
    (skip-chars-forward " \t")
    (let ((case-fold-search t))
      (if (or (looking-at "[0-9]")	;Reindent only where it is most
	      (looking-at "end")	;likely to be necessary
	      (looking-at "else")
	      (looking-at (regexp-quote fortran-continuation-string)))
	  (fortran-indent-line))))
  (newline)
  (fortran-indent-line))

(defun fortran-indent-subprogram ()
  "Properly indent the Fortran subprogram which contains point."
  (interactive)
  (save-excursion
    (mark-fortran-subprogram)
    (message "Indenting subprogram...")
    (indent-region (point) (mark) nil))
  (message "Indenting subprogram...done."))

(defun fortran-calculate-indent ()
  "Calculates the Fortran indent column based on previous lines."
  (let (icol first-statement (case-fold-search t)
	     (fortran-minimum-statement-indent
	      (if indent-tabs-mode
		  fortran-minimum-statement-indent-tab
		fortran-minimum-statement-indent-fixed)))
    (save-excursion
      (setq first-statement (fortran-previous-statement))
      (if first-statement
	  (setq icol fortran-minimum-statement-indent)
	(progn
	  (if (= (point) (point-min))
	      (setq icol fortran-minimum-statement-indent)
	    (setq icol (fortran-current-line-indentation)))
	  (skip-chars-forward " \t0-9")
	  (cond ((looking-at "\\(\\(\\sw\\|\\s_\\)+:[ \t]*\\)?if[ \t]*(")
		 (if (or (looking-at ".*)[ \t]*then\\b[ \t]*[^ \t_$(=a-z0-9]")
			 (let (then-test)	;multi-line if-then
			   (while (and (= (forward-line 1) 0)
				       ;;search forward for then
				       (or (looking-at "     [^ 0\n]")
					   (looking-at "\t[1-9]"))
				       (not (setq then-test (looking-at
							     ".*then\\b[ \t]\
*[^ \t_$(=a-z0-9]")))))
			   then-test))
		     (setq icol (+ icol fortran-if-indent))))
		((looking-at "else\\(if\\)?\\b")
		 (setq icol (+ icol fortran-if-indent)))
		((looking-at "select[ \t]*case[ \t](.*)")
		 (setq icol (+ icol fortran-if-indent)))
		((looking-at "case[ \t]*(.*)")
		 (setq icol (+ icol fortran-if-indent)))
		((looking-at "case[ \t]*default\\b")
		 (setq icol (+ icol fortran-if-indent)))
		((looking-at "\\(otherwise\\|else[ \t]*where\\)\\b")
		 (setq icol (+ icol fortran-if-indent)))
		((looking-at "where[ \t]*(.*)[ \t]*\n")
		 (setq icol (+ icol fortran-if-indent)))
		((looking-at "do\\b")
		 (setq icol (+ icol fortran-do-indent)))
		((looking-at
		  "\\(structure\\|union\\|map\\|interface\\)\\b[ \t]*[^ \t=(a-z]")
		 (setq icol (+ icol fortran-structure-indent)))
		((and (looking-at fortran-end-prog-re1)
		      (fortran-check-end-prog-re))
		 ;; Previous END resets indent to minimum
		 (setq icol fortran-minimum-statement-indent))))))
    (save-excursion
      (beginning-of-line)
      (cond ((looking-at "[ \t]*$"))
	    ((looking-at comment-line-start-skip)
	     (cond ((eq fortran-comment-indent-style 'relative)
		    (setq icol (+ icol fortran-comment-line-extra-indent)))
		   ((eq fortran-comment-indent-style 'fixed)
		    (setq icol (+ fortran-minimum-statement-indent
				  fortran-comment-line-extra-indent))))
	     (setq fortran-minimum-statement-indent 0))
	    ((or (looking-at (concat "[ \t]*"
				     (regexp-quote
				      fortran-continuation-string)))
		 (looking-at "     [^ 0\n]")
		 (looking-at "\t[1-9]"))
	     (setq icol (+ icol fortran-continuation-indent)))
	    ((looking-at "[ \t]*#")	; Check for cpp directive.
	     (setq fortran-minimum-statement-indent 0 icol 0))
	    (first-statement)
	    ((and fortran-check-all-num-for-matching-do
		  (looking-at "[ \t]*[0-9]+")
		  (fortran-check-for-matching-do))
	     (setq icol (- icol fortran-do-indent)))
	    (t
	     (skip-chars-forward " \t0-9")
	     (cond ((looking-at "end[ \t]*\\(if\\|select\\|where\\)\\b")
		    (setq icol (- icol fortran-if-indent)))
		   ((looking-at "else\\(if\\)?\\b")
		    (setq icol (- icol fortran-if-indent)))
                   ((looking-at "case[ \t]*\\((.*)\\|default\\>\\)")
		    (setq icol (- icol fortran-if-indent)))
		   ((looking-at "\\(otherwise\\|else[ \t]*where\\)\\b")
		    (setq icol (- icol fortran-if-indent)))
		   ((and (looking-at "continue\\b")
			 (fortran-check-for-matching-do))
		    (setq icol (- icol fortran-do-indent)))
		   ((looking-at "end[ \t]*do\\b")
		    (setq icol (- icol fortran-do-indent)))
		   ((looking-at "end[ \t]*\
\\(structure\\|union\\|map\\|interface\\)\\b[ \t]*[^ \t=(a-z]")
		    (setq icol (- icol fortran-structure-indent)))
		   ((and (looking-at fortran-end-prog-re1)
			 (fortran-check-end-prog-re)
			 (not (= icol fortran-minimum-statement-indent)))
 		    (message "Warning: `end' not in column %d.  Probably\
 an unclosed block." fortran-minimum-statement-indent))))))
    (max fortran-minimum-statement-indent icol)))

(defun fortran-current-line-indentation ()
  "Indentation of current line, ignoring Fortran line number or continuation.
This is the column position of the first non-whitespace character
aside from the line number and/or column 5/8 line-continuation character.
For comment lines, returns indentation of the first
non-indentation text within the comment."
  (save-excursion
    (beginning-of-line)
    (cond ((looking-at comment-line-start-skip)
	   (goto-char (match-end 0))
	   (skip-chars-forward
	    (if (stringp fortran-comment-indent-char)
		fortran-comment-indent-char
	      (char-to-string fortran-comment-indent-char))))
	  ((or (looking-at "     [^ 0\n]")
	       (looking-at "\t[1-9]"))
	   (goto-char (match-end 0)))
	  (t
	   ;; Move past line number.
	   (skip-chars-forward "[ \t0-9]");From Uli
	   ))
    ;; Move past whitespace.
    (skip-chars-forward " \t")
    (current-column)))

(defun fortran-indent-to-column (col)
  "Indent current line with spaces to column COL.
notes: 1) A non-zero/non-blank character in column 5 indicates a continuation
          line, and this continuation character is retained on indentation;
       2) If `fortran-continuation-string' is the first non-whitespace
          character, this is a continuation line;
       3) A non-continuation line which has a number as the first
          non-whitespace character is a numbered line.
       4) A TAB followed by a digit indicates a continuation line."
  (save-excursion
    (beginning-of-line)
    (if (looking-at comment-line-start-skip)
	(if fortran-comment-indent-style
	    (let ((char (if (stringp fortran-comment-indent-char)
			    (aref fortran-comment-indent-char 0)
			  fortran-comment-indent-char)))
	      (goto-char (match-end 0))
	      (delete-horizontal-regexp (concat " \t" (char-to-string char)))
	      (insert-char char (- col (current-column)))))
      (if (looking-at "\t[1-9]")
	  (if indent-tabs-mode
	      (goto-char (match-end 0))
	    (delete-char 2)
	    (insert "     ")
	    (insert fortran-continuation-string))
	(if (looking-at "     [^ 0\n]")
	    (if indent-tabs-mode
		(progn (delete-char 6)
		       (insert "\t")
		       (insert-char (fortran-numerical-continuation-char) 1))
	      (forward-char 6))
	  (delete-horizontal-space)
	  ;; Put line number in columns 0-4
	  ;; or put continuation character in column 5.
	  (cond ((eobp))
		((looking-at (regexp-quote fortran-continuation-string))
		 (if indent-tabs-mode
		     (progn
		       (indent-to
			(if indent-tabs-mode
			    fortran-minimum-statement-indent-tab
			  fortran-minimum-statement-indent-fixed))
		       (delete-char 1)
		       (insert-char (fortran-numerical-continuation-char) 1))
		   (indent-to 5)
		   (forward-char 1)))
		((looking-at "[0-9]+")
		 (let ((extra-space (- 5 (- (match-end 0) (point)))))
		   (if (< extra-space 0)
		       (message "Warning: line number exceeds 5-digit limit.")
		     (indent-to (min fortran-line-number-indent extra-space))))
		 (skip-chars-forward "0-9")))))
      ;; Point is now after any continuation character or line number.
      ;; Put body of statement where specified.
      (delete-horizontal-space)
      (indent-to col)
      ;; Indent any comment following code on the same line.
      (if (and comment-start-skip
	       (fortran-find-comment-start-skip))
	  (progn (goto-char (match-beginning 0))
		 (if (not (= (current-column) (fortran-comment-hook)))
		     (progn (delete-horizontal-space)
			    (indent-to (fortran-comment-hook)))))))))

(defun fortran-line-number-indented-correctly-p ()
  "Return t if current line's line number is correctly indented.
Do not call if there is no line number."
  (save-excursion
    (beginning-of-line)
    (skip-chars-forward " \t")
    (and (<= (current-column) fortran-line-number-indent)
	 (or (= (current-column) fortran-line-number-indent)
	     (progn (skip-chars-forward "0-9")
		    (= (current-column) 5))))))

(defun fortran-check-for-matching-do ()
  "When called from a numbered statement, return t if matching DO is found.
Otherwise return nil."
  (let (charnum
	(case-fold-search t))
    (save-excursion
      (beginning-of-line)
      (if (looking-at "[ \t]*[0-9]+")
	  (progn
	    (skip-chars-forward " \t")
	    (skip-chars-forward "0") ;skip past leading zeros
	    (setq charnum (buffer-substring (point)
					    (progn (skip-chars-forward "0-9")
						   (point))))
	    (beginning-of-line)
	    (fortran-with-subprogram-narrowing
	     (and (re-search-backward
		   (concat "\\(^[ \t0-9]*do[ \t]*0*" charnum "\\b\\)\\|"
			   "\\(^[ \t]*0*" charnum "\\b\\)")
		   nil t)
		  (looking-at (concat "^[ \t0-9]*do[ \t]*0*" charnum)))))))))

(defun fortran-find-comment-start-skip ()
  "Move to past `comment-start-skip' found on current line.
Return t if `comment-start-skip' found, nil if not."
  ;; In order to move point only if comment-start-skip is found, this
  ;; one uses a lot of save-excursions.  Note that re-search-forward
  ;; moves point even if comment-start-skip is inside a string-constant.
  ;; Some code expects certain values for match-beginning and end
  (interactive)
  (if (save-excursion
	(re-search-forward comment-start-skip
			   (save-excursion (end-of-line) (point)) t))
      (let ((save-match-beginning (match-beginning 0))
	    (save-match-end (match-end 0)))
	(if (fortran-is-in-string-p (match-beginning 0))
	    (save-excursion
	      (goto-char save-match-end)
	      (fortran-find-comment-start-skip)) ; recurse for rest of line
	  (goto-char save-match-beginning)
	  (re-search-forward comment-start-skip
			     (save-excursion (end-of-line) (point)) t)
	  (goto-char (match-end 0))
	  t))
    nil))

;;From: simon@gnu (Simon Marshall)
;; Find the next ! not in a string.
(defun fortran-match-!-comment (limit)
  (let (found)
    (while (and (setq found (search-forward "!" limit t))
                (fortran-is-in-string-p (point))))
    (if (not found)
	nil
      ;; Cheaper than `looking-at' "!.*".
      (set-match-data
       (list (1- (point)) (progn (end-of-line) (min (point) limit))))
      t)))

;; The above function is about 10% faster than the below...
;;(defun fortran-match-!-comment (limit)
;;  (let (found)
;;    (while (and (setq found (re-search-forward "!.*" limit t))
;;                (fortran-is-in-string-p (match-beginning 0))))
;;    found))

;;From: ralf@up3aud1.gwdg.de (Ralf Fassel)
;; Test if TAB format continuation lines work.
(defun fortran-is-in-string-p (where)
  "Return non-nil iff WHERE (a buffer position) is inside a Fortran string."
  (save-excursion
    (goto-char where)
    (cond
     ((bolp) nil)			; bol is never inside a string
     ((save-excursion			; comment lines too
	(beginning-of-line)
	(looking-at comment-line-start-skip)) nil)
     (t (let (;; ok, serious now. Init some local vars:
	      (parse-state '(0 nil nil nil nil nil 0))
	      (quoted-comment-start (if comment-start
					(regexp-quote comment-start)))
	      (not-done t)
	      parse-limit
	      end-of-line
	      )
	  ;; move to start of current statement
	  (fortran-next-statement)
	  (fortran-previous-statement)
	  ;; now parse up to WHERE
	  (while not-done
	    (if (or ;; skip to next line if:
		 ;; - comment line?
		 (looking-at comment-line-start-skip)
		 ;; - at end of line?
		 (eolp)
		 ;; - not in a string and after comment-start?
		 (and (not (nth 3 parse-state))
		      comment-start
		      (equal comment-start
			     (char-to-string (preceding-char)))))
		(if (> (forward-line) 0)
		    (setq not-done nil))
	      ;; else:
	      ;; if we are at beginning of code line, skip any
	      ;; whitespace, labels and tab continuation markers.
	      (if (bolp) (skip-chars-forward " \t0-9"))
	      ;; if we are in column <= 5 now, check for continuation char
	      (cond ((= 5 (current-column)) (forward-char 1))
		    ((and (< (current-column) 5)
			  (equal fortran-continuation-string
				 (char-to-string (following-char)))
			  (forward-char 1))))
	      ;; find out parse-limit from here
	      (setq end-of-line (save-excursion (end-of-line)(point)))
	      (setq parse-limit (min where end-of-line))
	      ;; parse max up to comment-start, if non-nil and in current line
	      (if comment-start
		  (save-excursion
		    (if (re-search-forward quoted-comment-start end-of-line t)
			(setq parse-limit (min (point) parse-limit)))))
	      ;; now parse if still in limits
	      (if (< (point) where)
		  (setq parse-state (parse-partial-sexp
				     (point) parse-limit nil nil parse-state))
		(setq not-done nil))
	      ))
	  ;; result is
	  (nth 3 parse-state))))))

(defun fortran-auto-fill-mode (arg)
  "Toggle fortran-auto-fill mode.
With ARG, turn `fortran-auto-fill' mode on iff ARG is positive.
In `fortran-auto-fill' mode, inserting a space at a column beyond `fill-column'
automatically breaks the line at a previous space."
  (interactive "P")
  (prog1 (setq auto-fill-function
	       (if (if (null arg)
		       (not auto-fill-function)
		     (> (prefix-numeric-value arg) 0))
		   #'fortran-do-auto-fill
		 nil))
    (force-mode-line-update)))

(defun fortran-do-auto-fill ()
  (if (> (current-column) fill-column)
      (fortran-indent-line)))

(defun fortran-fill ()
  (interactive)
  (let* ((auto-fill-function #'fortran-do-auto-fill)
	 (opoint (point))
	 (bol (save-excursion (beginning-of-line) (point)))
	 (eol (save-excursion (end-of-line) (point)))
	 (bos (min eol (+ bol (fortran-current-line-indentation))))
	 (quote
	  (save-excursion
	    (goto-char bol)
	    (if (looking-at comment-line-start-skip)
		nil			; OK to break quotes on comment lines.
	      (move-to-column fill-column)
	      (if (fortran-is-in-string-p (point))
		  (save-excursion (re-search-backward "\\S\"\\s\"\\S\"" bol t)
				  (if fortran-break-before-delimiters
				      (point)
				    (1+ (point))))))))
	 ;; decide where to split the line. If a position for a quoted
	 ;; string was found above then use that, else break the line
	 ;; before the last delimiter.
	 ;; Delimiters are whitespace, commas, and operators.
	 ;; Will break before a pair of *'s.
	 (fill-point
	  (or quote
	      (save-excursion
		(move-to-column (1+ fill-column))
		(skip-chars-backward "^ \t\n,'+-/*=)"
;;;		 (if fortran-break-before-delimiters
;;;		     "^ \t\n,'+-/*=" "^ \t\n,'+-/*=)")
		 )
		(if (<= (point) (1+ bos))
		    (progn
		      (move-to-column (1+ fill-column))
		      ;;what is this doing???
		      (if (not (re-search-forward "[\t\n,'+-/*)=]" eol t))
			  (goto-char bol))))
		(if (bolp)
		    (re-search-forward "[ \t]" opoint t)
		  (backward-char)
		  (if (looking-at "\\s\"")
		      (forward-char)
		    (skip-chars-backward " \t\*")))
		(if fortran-break-before-delimiters
		    (point)
		  (1+ (point)))))))
    ;; if we are in an in-line comment, don't break unless the
    ;; line of code is longer than it should be. Otherwise
    ;; break the line at the column computed above.
    ;;
    ;; Need to use fortran-find-comment-start-skip to make sure that quoted !'s
    ;; don't prevent a break.
    (if (not (or (save-excursion
		   (if (and (re-search-backward comment-start-skip bol t)
			    (not (fortran-is-in-string-p (point))))
		       (progn
			 (skip-chars-backward " \t")
			 (< (current-column) (1+ fill-column)))))
		 (save-excursion
		   (goto-char fill-point)
		   (bolp))))
	(if (> (save-excursion
		 (goto-char fill-point) (current-column))
	       (1+ fill-column))
	    (progn (goto-char fill-point)
		   (fortran-break-line))
	  (save-excursion
	    (if (> (save-excursion
		     (goto-char fill-point)
		     (current-column))
		   (+ (fortran-calculate-indent) fortran-continuation-indent))
		(progn
		  (goto-char fill-point)
		  (fortran-break-line))))))
    ))
(defun fortran-break-line ()
  (let ((opoint (point))
	(bol (save-excursion (beginning-of-line) (point)))
	(eol (save-excursion (end-of-line) (point)))
	(comment-string nil))

    (save-excursion
      (if (and comment-start-skip (fortran-find-comment-start-skip))
	  (progn
	    (re-search-backward comment-start-skip bol t)
	    (setq comment-string (buffer-substring (point) eol))
	    (delete-region (point) eol))))
    ;; Forward line 1 really needs to go to next non white line
    (if (save-excursion (forward-line)
			(or (looking-at "     [^ 0\n]")
			    (looking-at "\t[1-9]")))
	(progn
	  (end-of-line)
	  (delete-region (point) (match-end 0))
	  (delete-horizontal-space)
	  (fortran-fill))
      (fortran-split-line))
    (if comment-string
	(save-excursion
	  (goto-char bol)
	  (end-of-line)
	  (delete-horizontal-space)
	  (indent-to (fortran-comment-hook))
	  (insert comment-string)))))

(defun fortran-analyze-file-format ()
  "Return nil if fixed format is used, t if TAB formatting is used.
Use `fortran-tab-mode-default' if no non-comment statements are found in the
file before the end or the first `fortran-analyze-depth' lines."
  (let ((i 0))
    (save-excursion
      (goto-char (point-min))
      (setq i 0)
      (while (not (or
		   (eobp)
		   (looking-at "\t")
		   (looking-at "      ")
		   (> i fortran-analyze-depth)))
	(forward-line)
	(setq i (1+ i)))
      (cond
       ((looking-at "\t") t)
       ((looking-at "      ") nil)
       (fortran-tab-mode-default t)
       (t nil)))))

(or (assq 'fortran-tab-mode-string minor-mode-alist)
    (setq minor-mode-alist (cons
			    '(fortran-tab-mode-string
			      (indent-tabs-mode fortran-tab-mode-string))
			    minor-mode-alist)))

(defun fortran-fill-paragraph (&optional justify)
  "Fill surrounding comment block as paragraphs, else fill statement.

Intended as the value of `fill-paragraph-function'."
  (interactive "P")
  (save-excursion
    (beginning-of-line)
    (if (not (looking-at "[Cc*]"))
	(fortran-fill-statement)
      ;; We're in a comment block.  Find the start and end of a
      ;; paragraph, delimited either by non-comment lines or empty
      ;; comments.  (Get positions as markers, since the
      ;; `indent-region' below can shift the block's end).
      (let* ((non-empty-comment (concat "\\(" comment-line-start-skip
					"\\)" "[^ \t\n]"))
	     (start (save-excursion
		      ;; Find (start of) first line.
		      (while (and (zerop (forward-line -1))
				  (looking-at non-empty-comment)))
		      (or (looking-at non-empty-comment)
			  (forward-line)) ; overshot
		      (point-marker)))
	     (end (save-excursion
		    ;; Find start of first line past region to fill.
		    (while (progn (forward-line)
				  (looking-at non-empty-comment)))
		    (point-marker))))
	;; Indent the block, find the string comprising the effective
	;; comment start skip and use that as a fill-prefix for
	;; filling the region.
	(indent-region start end nil)
	(let ((paragraph-ignore-fill-prefix nil)
	      (fill-prefix (progn (beginning-of-line)
				  (looking-at comment-line-start-skip)
				  (match-string 0))))
	  (let (fill-paragraph-function)
	    (fill-region start end justify))) ; with normal `fill-paragraph'
	(set-marker start nil)
	(set-marker end nil))))
  t)

(defun fortran-fill-statement ()
  "Fill a fortran statement up to `fill-column'."
  (interactive)
  (let ((auto-fill-function #'fortran-do-auto-fill))
    (if (not (save-excursion
	       (beginning-of-line)
	       (or (looking-at "[ \t]*$")
		   (looking-at comment-line-start-skip)
		   (and comment-start-skip
			(looking-at (concat "[ \t]*" comment-start-skip))))))
	(save-excursion
	  ;; Find beginning of statement.
	  (fortran-next-statement)
	  (fortran-previous-statement)
	  ;; Re-indent initially.
	  (fortran-indent-line)
	  ;; Replace newline plus continuation field plus indentation with
	  ;; single space.
	  (while (progn
		   (forward-line)
		   (fortran-remove-continuation)))
	  (fortran-previous-statement)))
    (fortran-indent-line)))

; --- The following are extensions by Mark van Schilfgaarde ---
; 27 Jul 07 modified to handle limited f90 commands
; 14 Sep 07 New fortran-next-block and fortran-previous-block compatible with f90
; 23 Nov 11 New fortran-ifort-clean-unused-var spag-load-t-file fortran-to-std
(defun fortran-ifort-clean-unused-var()
" Uses output from Intel ifort to eliminate unused variables.  Run 'ifort -warn unused' in
 a shell and put the cursor at the beginning of the output."
  (interactive)

  (if (re-search-forward "This variable has not been used" (point-max) t)
      (let (vnam (fname (buffer-substring
			 (progn (beginning-of-line) (point))
			 (progn (re-search-forward "(") (- (point) 1))))
		 (lineno (buffer-substring
			  (point)
			  (progn (re-search-forward ")") (- (point) 1))))
		 )

	(find-file-other-window fname)
	(other-window -1)

	(while (re-search-forward "This variable has not been used" (point-max) t)
	  (re-search-forward "\\[\\([^]]+\\)")
	  (setq vnam (downcase (buffer-substring (match-beginning 1) (match-end 1))))
	  (message "name is %s" vnam)
	  (beginning-of-line)
;	  (setq lineno (progn (re-search-forward "(\\([0-9]+\\))") (match-string 1)))
	  (forward-line 1)
	  (other-window 1) (goto-line (string-to-int lineno)) (beginning-of-line)
	  (query-replace-regexp (concat ", *\\b" vnam "\\b\\(([^)]+)\\)*\\|\\b" vnam "\\b\\(([^)]+)\\)* *,\\|\\b" vnam "\\b\\(([^)]+)\\)*") "")
	  (other-window -1)))
    (message "string not found, 'This variable has not been used'")
    )
)


(defun spag-load-t-file ()
" Reads t.spg, puts in fortran mode, and rungs fortran-to-std"
  (interactive)

  (let ((bufnow (buffer-name))
	 strn)

    (delete-other-windows-quietly)
    (find-file-force "./t.spg" t)
    (split-window-quietly)
    (switch-to-buffer bufnow)
    (other-window 1)
    (fortran-mode)
    (goto-char (point-min))
    (fortran-to-std)
    )
)

(defun fortran-to-std ()
" Modifies spag output to my standard style"
  (interactive)

 (let  (decl-start                    ;not used
       (charvar "[a-z0-9_]+")         ;characters allowed in variable declarations
       (charexpr "[a-z0-9_*/()+-]+")  ;characters allowed in variable declarations
       )

   (setq decl-start (point))
;  (snit)

   (goto-char (point-min))
   (save-excursion
     (setq case-fold-search nil)
     (query-replace-regexp "\\([0-9]*[0-9.]+\\)\\([DE]\\)\\([+-]*[0-9]+\\)"  "\\1d\\3"))
   (goto-char (point-min))
   (query-replace-regexp "format +(" "format(")
   (query-replace-regexp "\\bif(" "if (")
   (goto-char (point-min))
   (query-replace-regexp (concat "\\([ \t]+do[ \t]\\)[ \t]*\\(" charvar "\\) *= *\\(" charexpr "\\) *, *\\(" charexpr "\\)") "\\1 \\2 = \\3, \\4")

   (goto-char (point-min))
   (while (re-search-forward " \\(continue\\|format *(\\)" (point-max) t)
     (progn
       (beginning-of-line)
       (if (looking-at "^ ") (indent-for-tab-command))
       (end-of-line)
       )
     )

   (goto-char (point-min))
   (query-replace-regexp "=\\([a-zA-Z0-9]\\)" " = \\1" nil nil nil)

   (goto-char decl-start)
))


(defun fortran-space-do-enddo ()
" Spaces out arguments of do-enddo construct"
  (interactive)
;  (snit)
;  (query-replace-regexp (concat "\\([ \t]+do[ \t]\\)[ \t]*\\(" charvar "\\) *= *\\(" charvar "\\) *, *\\(" charvar "\\)") "\\1 \\2 = \\3, \\4")

 (let  (decl-start                    ;not used
       (charvar "[a-z0-9_]+")         ;characters allowed in variable declarations
       (charexpr "[a-z0-9_*/()+-]+")  ;characters allowed in variable declarations
       )

  (query-replace-regexp (concat "\\([ \t]+do[ \t]\\)[ \t]*\\(" charvar "\\) *= *\\(" charexpr "\\) *, *\\(" charexpr "\\)") "\\1 \\2 = \\3, \\4")
))

(defvar fortran-start-of-prog-re1 "module\\|program\\|subroutine\\|[a-zA-Z]+[ \t]*function[ \t]*[a-zA-Z]+[ \t]*(")
(defvar fortran-end-of-prog-re1 (concat "end\\>[ \t]*" "\\(interface\\|function\\|module\\|program\\|subroutine\\)\\>" "\\|end\\>[ \t]*$"))
;(defvar fortran-do-start-re (concat "\\([ \t0-9]*do[ \t]*[0-9 \]*[A-za-z_ \t]+=\\)" "\\|\\([ \t0-9]*do[ \t]*while\\)"))
;(defvar fortran-do-start-re1 (concat "\\(do[ \t]*[0-9 \]*[A-za-z_ \t]+=\\)" "\\|\\(do[ \t]*while\\)"))

; --- macro fortran-first-executable-statement ---
;; Nov 2013 MvS updated to jump past "interface"
(defun fortran-first-executable-statement ()
"  Moves point to first executable statement in current subprogram."
   (interactive)

   (let ((case-fold-search t) pointnow (start (point)) )
   (setq subend (end-of-fortran-subprogram))
;   (snit)
   (beginning-of-fortran-subprogram)
   (re-search-forward "^[ 0-9][ 0-9][ 0-9][ 0-9][ 0-9][ 0-9] *") (beginning-of-line)
   (fortran-next-statement)
   (setq pointnow (point))
   (goto-char start)
   (f90-beginning-of-subprogram)
   (if (> pointnow (point)) (goto-char pointnow))

   ;; Nov 2013 handles INTERFACE as part of declaration
   (while (or (looking-at
            " +\\(character\\|common\\|use\\|type\\|complex\\|data\\|dimension\\|double\\|equivalence\\|external\\|function\\|implicit\\|integer\\|intrinsic\\|logical\\|parameter\\|program\\|real\\|save\\|subroutine\\|include\\)\\b")
	      (if (looking-at " +\\(interface\\)\\b") (re-search-forward "^[ \t]+end[ \t]+interface" subend nil)))
      (fortran-next-statement))
))

(defun get-query (msg unit)
" Returns one of y, n or a, and queries if unit not a"
(interactive)
  (if (not (char-equal unit ?a))
     (progn (message msg) (setq unit (read-char)))
     (setq unit unit))
)

(defun fortran-const-to-double (startposn posn que)
  "converts fortran constants to double precision in region (STARTPOSN, POSN).
   If QUE (third arg) is t, asks interactively.
   Returns change in size of region."
(interactive)
; Set up for query replace
    (if que (setq query ??) (setq query ?a))
    (setq savpos posn)
    (goto-char startposn)
;   The following is a necessary condition for string to be a real number
    (while (re-search-forward "[0-9]" posn t)
       (if
; Exclude strings beginning with [a-z]
	     (save-excursion
	       (backward-char 1) (skip-chars-backward "A-Za-z ")
               (while (< (current-column) 7) (forward-char 1))
	       (not (looking-at " *[a-zA-Z]")))
; OK, we have a bonafide number
         (let ((case-fold-search t))
; Determine whether number contains a "."
	   (setq string_contains_point
	     (or
	       (and (looking-at "[0-9 ]*[.]")
                    (not (looking-at
			   "[0-9 ]*[.]\\(eq\\|and\\|or\\|gt\\|ge\\|lt\\|le\\|ne\\)[.]")))
	     (save-excursion
	       (backward-char 1)
	       (skip-chars-backward "A-Za-z .")
               (while (< (current-column) 7) (forward-char 1))
	       (skip-chars-forward "A-Za-z ")
	       (and (looking-at " *[.]")
                    (not (looking-at
			   " *[.]\\(eq\\|and\\|or\\|gt\\|ge\\|lt\\|le\\|ne\\)[.]"))))))
; Skip past mantissa
	   (skip-chars-forward "0-9 ") (skip-chars-forward ".") (skip-chars-forward "0-9 ")
; If already double-precision, merely skip past exponent, otherwise  ...
          (if (looking-at "[Dd]") (skip-chars-forward "-+0-9") (progn
; Convert real with explicit exponent: string ensures exponent a pure integer
             (if (looking-at "[Ee][-+0-9]+\\b[^.]")
	        (progn
		   (setq query (get-query "Convert?" query))
		   (if (or (char-equal query ?a) (char-equal query ?y))
; either replace it or skip over it ...
		       (progn (delete-char 1) (insert-char ?d 1)) (forward-char 1))
; and then skip past exponent ...
		   (skip-chars-forward "-+0-9"))
; Convert real without explicit exponent
		(if string_contains_point (progn
		   (skip-chars-backward " ")
		   (setq query (get-query "Convert?" query))
		   (if (or (char-equal query ?a) (char-equal query ?y))
		       (progn (insert-string "d0") (setq posn (+ posn 2))))
      ))))))
; If not a number, skip past digits to avoid confusion about subsequent while
      (skip-chars-forward "0-9. ")))
(- posn savpos))                    ; Return with change in region size

(defun fortran-decl-to-double (startposn posn que)
  "converts fortran cast declarations to dble in region (STARTPOSN, POSN).
   If QUE (third arg) is t, asks interactively.
   Returns change in size of region."
(interactive)
; Set up for query replace
   (if que (setq query ??) (setq query ?a))
   (setq savpos posn)

; Handle case REAL is a declaration
   (goto-char startposn)
   (while (re-search-forward "^ +real\\b\\(\\*[0-9]+\\)* +[^(]" posn t) (progn
      (setq query (get-query "Convert?" query))
      (if (or (char-equal query ?a) (char-equal query ?y)) (progn
;        This kills the REAL
         (beginning-of-line) (forward-word 1) (backward-kill-word 1)
;        This kills the *number declaration
         (while (looking-at "[*0-9 ]")
            (progn (delete-char 1) (setq posn (- posn 1))))
            (insert-string "double precision")
            (insert-char ?\  1) (setq posn (+ posn 1))
            (setq posn (+ posn 12))
   ))))

; Handle case COMPLEX is a declaration
   (goto-char startposn)
   (while (re-search-forward "^ +complex\\b\\(\\*[0-9]+\\)* +[^(]" posn t) (progn
      (setq query (get-query "Convert?" query))
      (if (or (char-equal query ?a) (char-equal query ?y)) (progn
;        This kills the COMPLEX
         (beginning-of-line) (forward-word 1) (backward-kill-word 1)
;        This kills the *number declaration
         (while (looking-at "[*0-9 ]") (progn (delete-char 1) (setq posn (- posn 1))))
         (insert-string "complex*16")
         (insert-char ?\  1) (setq posn (+ posn 1))
         (setq posn (+ posn 3))
   ))))

; Handle implicit REAL declaration
   (goto-char startposn)
   (while (re-search-forward "^ +implicit +real\\b\\(\\*[0-9]+\\)* +" posn t) (progn
      (setq query (get-query "Convert?" query))
      (if (or (char-equal query ?a) (char-equal query ?y)) (progn
;        This kills the REAL
         (beginning-of-line) (forward-word 2) (backward-kill-word 1)
;        This kills the *number declaration
         (while (looking-at "[*0-9 ]")
         (progn (delete-char 1) (setq posn (- posn 1))))
	 (insert-string "double precision")
	 (insert-char ?\  1) (setq posn (+ posn 1))
	 (setq posn (+ posn 12))
   ))))

; Handle implicit COMPLEX declaration
   (goto-char startposn)
   (while (re-search-forward "^ +implicit +complex\\b\\(\\*[0-9]+\\)* +" posn t) (progn
      (setq query (get-query "Convert?" query))
      (if (or (char-equal query ?a) (char-equal query ?y)) (progn
;     This kills the COMPLEX
      (beginning-of-line) (forward-word 2) (backward-kill-word 1)
;     This kills the *number declaration
      (while (looking-at "[*0-9 ]")
         (progn (delete-char 1) (setq posn (- posn 1))))
         (insert-string "complex*16")
         (insert-char ?\  1) (setq posn (+ posn 1))
         (setq posn (+ posn 3))
   ))))

 (- posn savpos))                    ; Return with change in region size

(defun fortran-intrsc-to-double (startposn posn que)
  "converts intrinsic fortran functions to dble in region (STARTPOSN, POSN).
   If QUE (third arg) is t, asks interactively.
   Returns change in size of region."
  (interactive)
; Set up for query replace
   (if que (setq query ??) (setq query ?a))
   (setq savpos posn)
; Handle all cases where fnct -> dfnct
   (goto-char startposn)
   (while (re-search-forward
           "\\b\\(abs\\|acos\\|asin\\|atan\\|atan2\\|cmplx\\|conjg\\|cos\\|cosh\\|exp\\|float\\|sign\\|sin\\|sinh\\|sqrt\\|tan\\|tanh\\)\\b"  posn t) (progn
      (setq query (get-query "Convert?" query))
      (if (or (char-equal query ?a) (char-equal query ?y))
         (progn (backward-word 1) (insert-string "d")
            (setq posn (+ posn 1))
   ))))

; Handle all cases where cfnct -> cdfnct
    (goto-char startposn)
    (while (re-search-forward
            "\\bc\\(abs\\|cos\\|exp\\|log\\|sin\\|tan\\)\\b"  posn t) (progn
      (setq query (get-query "Convert?" query))
      (if (or (char-equal query ?a) (char-equal query ?y))
         (progn (backward-word 1) (forward-char 1) (insert-string "d")
              (setq posn (+ posn 1))
   ))))

; Handle all cases where afnct -> dfnct
   (goto-char startposn)
   (while (re-search-forward
           "\\ba\\(imag\\|log10\\|log\\|max1\\|min1\\|mod\\)\\b" posn t) (progn
      (setq query (get-query "Convert?" query))
      (if (or (char-equal query ?a) (char-equal query ?y))
         (progn (backward-word 1) (delete-char 1) (insert-string "d"))
   )))

; Handle SNGL -> DBLE
   (goto-char startposn)
   (while (re-search-forward "\\bsngl\\b" posn t) (progn
      (setq query (get-query "Convert?" query))
      (if (or (char-equal query ?a) (char-equal query ?y))
         (progn (backward-kill-word 1) (insert-string "dble"))
   )))

; Handle REAL -> DREAL
   (goto-char startposn)
   (while (re-search-forward "\\breal\\b" posn t)
      (if (looking-at " *(") (progn
      (setq query (get-query "Convert?" query))
      (if (or (char-equal query ?a) (char-equal query ?y))
         (progn (backward-word 1) (insert-string "d") (setq posn (+ posn 1)))
   ))))
 (- posn savpos))                    ; Return with change in region size

(defun fortran-to-double (startposn posn &optional que)
  "converts constants, declarations and intrinsic functions
   to double precision in region (STARTPOSN, POSN).
   If QUE (third arg) is t, asks interactively."
(interactive)
  (if (< startposn posn) (progn
    (setq posn (+ posn (fortran-const-to-double startposn posn que)))
    (setq posn (+ posn (fortran-decl-to-double startposn posn que)))
    (setq posn (+ posn (fortran-intrsc-to-double startposn posn que)))
    (goto-char posn)
)))

(defun fortran-cnvt-buffer-to-double ()
  "Converts constants, declarations and intrinsic functions to double
   precision in current buffer"
(interactive)
   (save-excursion
     (untabify (point-min) (point-max))
     (beginning-of-buffer)
     (while (< (point) (point-max))
; skip over comments
        (while (looking-at "^[Cc*]") (next-line 1))
        (fortran-to-double (point) (progn (end-of-line) (point)) t)
        (forward-char 1)
     )
))

(defun fortran-intrsc-to-sngl (startposn posn que)
  "converts intrinsic fortran functions to sngl in region (STARTPOSN, POSN).
   If QUE (third arg) is t, asks interactively.
   Returns change in size of region."
  (interactive)
; Set up for query replace
   (if que (setq query ??) (setq query ?a))
   (setq savpos posn)
; Handle all cases where dfnct -> fnct
   (goto-char startposn)
   (while (re-search-forward
           "\\bd\\(abs\\|acos\\|asin\\|atan\\|atan2\\|cmplx\\|conjg\\|cos\\|cosh\\|exp\\|float\\|real\\|sign\\|sin\\|sinh\\|sqrt\\|tan\\|tanh\\)\\b"  posn t) (progn
      (setq query (get-query "Convert?" query))
      (if (or (char-equal query ?a) (char-equal query ?y))
         (progn (backward-word 1) (delete-char 1)
            (setq posn (- posn 1))
   ))))

; Handle all cases where cdfnct -> cfnct
    (goto-char startposn)
    (while (re-search-forward
            "\\bcd\\(abs\\|cos\\|exp\\|log\\|sin\\|tan\\)\\b"  posn t) (progn
      (setq query (get-query "Convert?" query))
      (if (or (char-equal query ?a) (char-equal query ?y))
         (progn (backward-word 1) (forward-char 1) (delete-char 1)
              (setq posn (- posn 1))
   ))))

; Handle all cases where dfnct -> afnct
   (goto-char startposn)
   (while (re-search-forward
           "\\bd\\(imag\\|log10\\|log\\|max1\\|min1\\|mod\\)\\b" posn t) (progn
      (setq query (get-query "Convert?" query))
      (if (or (char-equal query ?a) (char-equal query ?y))
         (progn (backward-word 1) (delete-char 1) (insert-string "a"))
   )))

; Handle DBLE -> SNGL
   (goto-char startposn)
   (while (re-search-forward "\\bdble\\b" posn t) (progn
      (setq query (get-query "Convert?" query))
      (if (or (char-equal query ?a) (char-equal query ?y))
         (progn (backward-kill-word 1) (insert-string "sngl"))
   )))

 (- posn savpos))                    ; Return with change in region size


; --- macro fortran-compare-call-args ---
(defun fortran-compare-call-args ()
"   Compare arguments of a subroutine call with subroutine args"
   (interactive)

; ... Setup
   (let (thisbuf subnam substrng strngsiz thispt endpt subrstrng
	 sizcol r start end)
     (beginning-of-line)
     (fortran-next-statement)
     (fortran-previous-statement)
     (setq thisbuf (current-buffer))
     (setq tagnam (find-tag-default-fortran-mode))
     (setq subnam (concat tagnam "\\b"))
     (re-search-forward subnam nil 0)

; ... list all args of sub
   (setq sizcol 0)
   (save-excursion
     (find-tag-regexp subnam)
     (re-search-forward subnam nil 0)
     (skip-chars-forward " \t")
     (setq subrstrng (buffer-substring (+ 1 (point)) (progn (forward-sexp) (- (point) 1))))
     (pop-to-buffer "*compare*") (erase-buffer) (insert subrstrng "\n")
; .. Clean up around newlines
     (goto-char (point-min))
     (while (progn (forward-line 1) (< (point) (point-max)))
	(delete-char 6)
	(while (looking-at "[ \t]") (delete-char 1))
	(insert ". "))
; .. Put arguments into columns
     (goto-char (point-min))
     (while (< (point) (point-max))
       (forward-sexp) (setq sizcol (max sizcol (current-column)))
       (while (looking-at "[ \n\t,]") (delete-char 1)) (insert "\n")))
   (setq sizcol (+ sizcol 4))

; ... list all args of call
   (setq subrstrng (buffer-substring
		     (+ 1 (point))
		     (progn (forward-sexp) (- (point) 1))))
   (setq strngsiz 1)
   (pop-to-buffer "*compare*") (goto-char (point-min))
   (while (> strngsiz 0)
; ... Replace "^ \." with a newline now
     (beginning-of-line)
     (if (looking-at "^\\.") (progn (delete-char 2) (insert "\n")))
     (end-of-line) (while (< (current-column) sizcol) (insert " "))
     (setq thispt (point)) (insert subrstrng) (setq endpt (point))
     (goto-char thispt) (fortran-forward-expr 1 3)
; ... Skip past argument separators
     (while (and (looking-at "[ \t,]") (> endpt (point)))
       (delete-char 1) (setq endpt (- endpt 1)))
; ... Clear garbage around fortran continuation line
     (if (and (looking-at "\n") (> endpt (point)))
	 (progn (delete-char 7) (setq endpt (- endpt 7))
		(while (and (looking-at "[ \t,]") (> endpt (point)))
		  (delete-char 1) (setq endpt (- endpt 1)))))
     (setq subrstrng (buffer-substring (point) endpt))
     (setq strngsiz (- endpt (point)))
     (delete-char (- endpt (point)))
     (next-line 1))

; ... Mark all arguments with different names
   (pop-to-buffer "*compare*") (goto-char (point-min))
   (while (< (point) (point-max))
     (if (not (or (looking-at "^\\(.+\\) +\\1$")
		  (looking-at "^$")))
	 (progn (end-of-line) (insert " *")))
     (end-of-line) (setq sizcol (max sizcol (current-column)))
     (forward-line))

   (setq sizcol (+ sizcol 8))

; ... Double up columns if needed and room available
   (if (and
	 (>= (count-lines 1 (point-max)) (window-height))
	 (>= (window-width) (+ sizcol sizcol)))
       (progn
	 (goto-char (point-min)) (end-of-line)
	 (while (< (current-column) (/ (window-width) 2)) (insert " "))
	 (setq start
	       (progn (goto-line
			(+ (/ (+ (count-lines 1 (point-max)) 1) 2) 1))
		      (point)))
	 (setq end
	       (progn (goto-char (point-max))
		      (forward-line -1) (end-of-line)
		      (while (< (current-column) sizcol) (insert " "))
		      (point)))
	 (goto-char (point-min)) (end-of-line)
	 (insert-rectangle (delete-extract-rectangle start end))
	 ))

   (goto-char (point-min))
   )
)

; --- macro fortran-previous-block ---
(defun fortran-previous-block (&optional arg)
" If current-line is at the END of a DO- loop, go to of beginning-of-loop
 If current-line is at an ENDIF or ELSE or ELSEIF statement,
   go to to preceding ELSE or ELSE-IF it exists; otherwise go to the matching IF-THEN
 Otherwise, execute fortran-previous-statement until a beginning-of-block or end-of-block is found.
 Optional C-U causes this last option to be executed regardless of what line point is initially on.
  When invoked from start-of-block, causes point to move down a nesting level."

    (interactive "P")

    (let ((ifboundary (fortran-check-and-find-block-boundary 3))
	  elsepos narg priorif priordo pointnow
	  (prefix-is-C-u (and arg (not (equal arg (prefix-numeric-value arg))))))

      ; Set narg to :
      ; -1 if C-u
      ; -2 if C-u C-u
      ; -4 if C-u C-u C-u
      ;  0 if arg is nil
      ; arg if arg is a numerical value
      (if prefix-is-C-u
	  (if (= 4 (prefix-numeric-value arg))
	      (setq narg -1)
	    (if (= 16 (prefix-numeric-value arg))
		(setq narg -2)
	      (setq narg -4)))
        (if arg (setq narg (prefix-numeric-value arg))
	  (setq narg 0)))

      (if (or (not ifboundary)  (= narg -1))
	  ;; We are not at a start-boundary; move up until we find either start- or end- boundary
	  (while (not
		  (or (fortran-previous-statement) (fortran-check-and-find-block-boundary 1))))

        ;; We are at a start-boundary; move to corresponding end-boudary
	(progn

	  ;; If boundary is the start of if-then block, shift end to matching else or endif, if closer
	  (beginning-of-line) (skip-chars-forward " \t0-9")
	  (if (and

	       (or
		(looking-at "else[ 	]*if\\b")
		(looking-at "else\\b")
		(looking-at "end[ \t]*if\\b"))

	       ;; find possible following matching ELSE that is closer than IF-THEN
	       (setq elsepos (fortran-beginning-if-or-else)))
	      (setq ifboundary (max ifboundary elsepos)))

	  (goto-char ifboundary)

	  ;; Handle C-U C-U
	  (if (= narg -2)
	      (cond
	       ;; case looking at  at "else": or "elseif"
	       ((or
		   (looking-at "else[ 	]*if\\b")
		   (looking-at "else\\b"))
		     (goto-char (fortran-beginning-if)))

	       ((or (save-excursion (fortran-previous-statement) (fortran-check-and-find-block-boundary 3)))
		(progn

		  ;; move up any blocks at same nest level as starting point
		  (while (progn (fortran-previous-statement) (fortran-check-and-find-block-boundary 3))
		    (goto-char (fortran-check-and-find-block-boundary 3)))

		  ;; move up one nesting level, if present

					;set priordo to max(priorif,priordo) if both exist, otherwise to the value of whichever exists (nil if neither)
		  (setq priordo (fortran-beginning-do))
		  (setq priorif (fortran-beginning-if))
;		  (snit)

		  (if (not priordo) (setq priordo priorif))
		  (and priorif priordo (setq priordo (max priorif priordo)))

					;goto to bottom of prior block, if above current point; otherwise to top of block
		  (if priordo
		      (progn
			(setq pointnow (point))
			(goto-char priordo)
			(if (and (fortran-check-and-find-block-boundary 2) (< (fortran-check-and-find-block-boundary 2) pointnow))
			    (goto-char (fortran-check-and-find-block-boundary 2))))

		    (fortran-next-statement)) ))))

	    (beginning-of-line)))
))

; --- macro fortran-next-block ---
(defun fortran-next-block (&optional arg)
" If current-line is at the BEGINNING of a DO- loop, go to of end-of-loop
 If current-line is at an IF-THEN or ELSE or ELSEIF statement,
   go to to next ELSE or ELSE-IF it exists; otherwise go to the matching ENDIF.
 Otherwise, execute fortran-next-statement until a beginning-of-block or end-of-block is found.
 Optional C-U causes this last option to be executed regardless of what line point is initially on.
  When invoked from start-of-block, causes point to move down a nesting level."

    (interactive "P")

    (let ((ifboundary (fortran-check-and-find-block-boundary 2))
	  elsepos narg priorif priordo pointnow
	  (prefix-is-C-u (and arg (not (equal arg (prefix-numeric-value arg))))))

      ; Set narg to :
      ; -1 if C-u
      ; -2 if C-u C-u
      ; -4 if C-u C-u C-u
      ;  0 if arg is nil
      ; arg if arg is a numerical value
      (if prefix-is-C-u
	  (if (= 4 (prefix-numeric-value arg))
	      (setq narg -1)
	    (if (= 16 (prefix-numeric-value arg))
		(setq narg -2)
	      (setq narg -4)))
        (if arg (setq narg (prefix-numeric-value arg))
	  (setq narg 0)))

      (if (or (not ifboundary)  (= narg -1))
	  ;; We are not at a start-boundary; move down until we find either start- or end- boundary
	  (while (not
		  (or (fortran-next-statement) (fortran-check-and-find-block-boundary 1))))

        ;; We are at a start-boundary; move to corresponding end-boudary
	(progn

	  ;; If boundary is the start of if-then block, shift end to matching else or endif, if closer
	  (beginning-of-line) (skip-chars-forward " \t0-9")
	  (if (and

	       (or
		(looking-at "else[ 	]*if\\b")
		(looking-at "else\\b")
		(looking-at fortran-if-start-re))

	       ;; find possible following matching ELSE that is closer than ENDIF
	       (setq elsepos (fortran-end-if-or-else 2)))
	      (setq ifboundary (min ifboundary elsepos)))

	  (goto-char ifboundary)

	  ;; Handle C-U C-U
	  (if (= narg -2)
	      (cond
	       ;; case looking at  at "else" or "elseif"
	       ((or
		   (looking-at "else[ 	]*if\\b")
		   (looking-at "else\\b"))
		     (goto-char (fortran-end-if)))

	       ((or (save-excursion (fortran-next-statement) (fortran-check-and-find-block-boundary 2)))
		(progn

		  ;; move down any continguous blocks at same nest level as starting point
		  (while (progn (fortran-next-statement) (fortran-check-and-find-block-boundary 2))
		    (goto-char (fortran-check-and-find-block-boundary 2)))

		  ;; move up one nesting level, if present

					;set priordo to max(priorif,priordo) if both exist, otherwise to the value of whichever exists (nil if neither)
		  (setq priordo (fortran-end-do))
		  (setq priorif (fortran-end-if))

		  (if (not priordo) (setq priordo priorif))
		  (and priorif priordo (setq priordo (min priorif priordo)))

					;goto to top of prior block, if below current point; otherwise to bottomn of block
		  (if priordo
		      (progn
			(setq pointnow (point))
			(goto-char priordo)
			(if (and (fortran-check-and-find-block-boundary 3) (> (fortran-check-and-find-block-boundary 3) pointnow))
			    (goto-char (fortran-check-and-find-block-boundary 3))))

		    (fortran-previous-statement)) ))))

	    (beginning-of-line)))
))

; --- macro fortran-check-and-find-block-boundary ---
(defun fortran-check-and-find-block-boundary (&optional arg)
" If current-line is at the BEGINNING of either a DO- or IF-THEN block or subroutine, return of end-block.
 If current-line is at the END of either a DO- or IF-THEN block or subroutine, return of beginning-of-block.
 Optional ARG:
   If 2, include only cases for starting blocks: IF-THEN and DO
   If 3, include only cases for ending blocks: ENDIF or end of do loop."
   (interactive "p")


   (let ((start-ok (not (= arg 3)))
	 (end-ok (not (= arg 2))))

       (save-excursion
	 (beginning-of-line) (skip-chars-forward " \t0-9")
	 (cond

	  ;; start or middle of IF-THEN block: move to end of it
	  ((and start-ok
		(or
					; and block checks for if-then, including multi-line case
	    (and
	     (looking-at fortran-if-start-re)
	     (save-match-data
	       (or (looking-at ".*)[ \t]*then\\b[ \t]*[^ \t(=a-z0-9]")
		   ;; Multi-line if-then.
		   (save-excursion
		     (let (then-test)
		       (while
			   (and (= (forward-line 1) 0)
				;; Search forward for then.
				(or (looking-at "     [^ 0\n]")
				    (looking-at "\t[1-9]"))
				(not
				 (setq then-test
				       (looking-at
					".*then\\b[ \t]*[^ \t(=a-z0-9]")))))
		       then-test)))))
	    (looking-at "else[ 	]*if\\b")
	    (looking-at "else\\b")))
	   (fortran-end-if))

	  ;; end of IF-THEN block: move to start of it
	  ((and end-ok
		(or (looking-at "end[ \t]*if\\b")
		    (looking-at "else[ 	]*if\\b")
		    (looking-at "else\\b")))
	   (fortran-beginning-if))

	  ;; start of do-loop: move to end of it
	  ((and start-ok (fortran-check-for-matching-enddo))
	   (save-excursion (goto-char (match-beginning 0)) (skip-chars-forward " \t0-9") (point)))

	  ;; end of do-loop: move to start of it
	  ((and end-ok (fortran-check-for-matching-do))
	   (save-excursion (goto-char (match-beginning 0)) (skip-chars-forward " \t0-9") (point)))

	  ;; start-of-program
	  ((and start-ok (looking-at fortran-start-of-prog-re1))
	   (save-excursion (end-of-fortran-subprogram) (backward-char 1) (beginning-of-line) (point)))

	  ;; end-of-program
	  ((and end-ok (looking-at fortran-end-of-prog-re1))
	   (save-excursion (beginning-of-fortran-subprogram) (point)))

))))

; --- macro fortran-check-for-block-boundary-word.  Not particularly useful ---
; (defun fortran-check-for-block-boundary-word ()
; "If current-line is sitting on a keyword denoting beginning or end of either a DO- or IF-THEN block,
;  then return point.  If not, return nil."
;    (interactive)
;    (if
;        (save-excursion (beginning-of-line)
; 		       (skip-chars-forward " \t0-9")
; 		       (or
; 			(looking-at "end[ \t]*if\\b")
; 			(looking-at "else[ 	]*if\\b")
; 			(looking-at "else\\b")
; 			;; check for if-then, single or multiline case
; 			(and
; 			 (looking-at fortran-if-start-re)
; 			 (save-match-data
; 			   (or (looking-at ".*)[ \t]*then\\b[ \t]*[^ \t(=a-z0-9]")
; 			       ;; Multi-line if-then.
; 			       (save-excursion
; 				 (let (then-test)
; 				   (while
; 				       (and (= (forward-line 1) 0)
; 					    ;; Search forward for then.
; 					    (or (looking-at "     [^ 0\n]")
; 						(looking-at "\t[1-9]"))
; 					    (not
; 					     (setq then-test
; 						   (looking-at
; 						    ".*then\\b[ \t]*[^ \t(=a-z0-9]")))))
; 				   then-test)))))
; 			(if (fortran-check-for-matching-enddo) (looking-at "."))
; 			(if (fortran-check-for-matching-do) (looking-at "."))
; ; 			(looking-at fortran-do-start-re1)
; 			(looking-at fortran-end-prog-re1)))
;        ;; Sitting on one.
;        (match-beginning 0)
;      nil)
; )

; --- macro fortran-check-for-matching-enddo ---
(defun fortran-check-for-matching-enddo ()
  "When called from a numbered do statement, returns point of loop end if found, otherwise nil.
   When called from an numbered do statement, returns point of matching enddo if found, otherwise nil."
  (let (charnum
	(case-fold-search t))

    (save-excursion
      (beginning-of-line)
      (cond

       ;; numbered do
       ((looking-at "[ \t0-9]*do[ \t]*\\([0-9]+\\)[ \t]*[A-za-z_ \t]+=")
	(progn
	  (goto-char (match-beginning 1))
	  (setq charnum (buffer-substring (point)
					  (progn (skip-chars-forward "0-9")
						 (point))))
	  (re-search-forward (concat "^[ \t]*0*" charnum "\\b") nil t)))

       ;; do-enddo
      ((looking-at (concat "\\([ \t0-9]*do[ \t]*[A-za-z0-9_ \t]+=\\)" "\\|\\([ \t0-9]*do[ \t]*while\\)"))
       (fortran-end-do))
)
)))

; --- macro fortran-check-for-matching-do ---
(defun fortran-check-for-matching-do ()
  "When called from a numbered statement, returns point of matching do if found, otherwise nil.
   When called from an enddo statement, returns point of matching do if found, otherwise nil."
  (let (charnum
	(case-fold-search t))

    (save-excursion
      (beginning-of-line)
      (cond

       ;; numbered do
       ((looking-at "[ \t]*[0-9]+")
	(progn
	  (skip-chars-forward " \t")
	  (skip-chars-forward "0") ;skip past leading zeros
	  (setq charnum (buffer-substring (point)
					  (progn (skip-chars-forward "0-9")
						 (point))))
	  (beginning-of-line)
	  (and (re-search-backward
		(concat "\\(^[ \t0-9]*end\\b[ \t]*[^ \t=(a-z]\\)\\|"
			"\\(^[ \t0-9]*do[ \t]*0*" charnum "\\b\\)\\|"
			"\\(^[ \t]*0*" charnum "\\b\\)")
		nil t)
	       (if (looking-at (concat "^[ \t0-9]*do[ \t]*0*" charnum)) (match-beginning 0) nil))))

       ((looking-at "^[ \t0-9]+end[ \t]*do\\b")
	(fortran-beginning-do))

       ))))

; --- macro fortran-check-for-matching-numbered-do ---
(defun fortran-check-for-matching-numbered-do ()
  "When called from a numbered statement, Returns point if found, otherwise nil."
  (let (charnum
	(case-fold-search t))

    (save-excursion
      (beginning-of-line)
      (if (looking-at "[ \t]*[0-9]+")
	    (progn
	      (skip-chars-forward " \t")
	      (skip-chars-forward "0") ;skip past leading zeros
	      (setq charnum (buffer-substring (point)
					      (progn (skip-chars-forward "0-9")
						     (point))))
	      (beginning-of-line)
	      (and (re-search-backward
		    (concat "\\(^[ \t0-9]*end\\b[ \t]*[^ \t=(a-z]\\)\\|"
			    "\\(^[ \t0-9]*do[ \t]*0*" charnum "\\b\\)\\|"
			    "\\(^[ \t]*0*" charnum "\\b\\)")
		    nil t)
		   (if (looking-at (concat "^[ \t0-9]*do[ \t]*0*" charnum)) (match-beginning 0) nil)))
))))

; --- macro fortran-beginning-if-or-else ---
(defun fortran-beginning-if-or-else ()
  "Search backwards for first unmatched IF-THEN or ELSE-THEN. Return point or nil."
; Optional ARG=2 => if sitting at an ELSE, do not count it."

 (interactive)

  (let ((case-fold-search t))
    (if (save-excursion
	  ;; May be sitting on multi-line if-then statement, first move to
	  ;; beginning of current statement.  Note: `fortran-previous-statement'
	  ;; moves to previous statement *unless* current statement is first
	  ;; one.  Only move forward if not first-statement.
	  (if (not (eq (fortran-previous-statement) 'first-statement))
	      (fortran-next-statement))
	  (skip-chars-forward " \t0-9")
	  (and
	   (looking-at fortran-if-start-re)
	   (save-match-data
	     (or (looking-at ".*)[ \t]*then\\b[ \t]*[^ \t(=a-z0-9]")
		 ;; Multi-line if-then.
		 (let (then-test)
		   (while
                     (and (= (forward-line 1) 0)
			    ;; Search forward for then.
			    (or (looking-at "     [^ 0\n]")
				(looking-at "\t[1-9]"))
			    (not
			     (setq then-test
				   (looking-at
				    ".*then\\b[ \t]*[^ \t(=a-z0-9]")))))
		   then-test)))))
	;; Sitting on one.
	(match-beginning 0)
      ;; Search for one.
      (save-excursion
	(let ((count 1))
        (while (and (not (= count 0))
		      (not (eq (fortran-previous-statement) 'first-statement))
		      ;; Keep local to subprogram.
		      (not (and (looking-at fortran-end-prog-re)
				(fortran-check-end-prog-re))))

	    (skip-chars-forward " \t0-9")
	    (cond
	     ((looking-at fortran-if-start-re)
		   (save-excursion
		     (if (or
			  (looking-at ".*)[ \t]*then\\b[ \t]*[^ \t(=a-z0-9]")
			  (let (then-test) ; Multi-line if-then.
			    (while
                              (and (= (forward-line 1) 0)
				     ;; Search forward for then.
				     (or (looking-at "     [^ 0\n]")
					 (looking-at "\t[1-9]"))
				     (not
				      (setq then-test
					    (looking-at
					     ".*then\\b[ \t]*[^ \t(=a-z0-9]")))))
			    then-test))
                       (setq count (- count 1)))))

	     ((and (= count 1) (looking-at "\\(else[ 	]*if\\b\\)\\|\\(else\\b\\)"))
	      (setq count (- count 1)))

	     ((looking-at "end[ \t]*if\\b")
	      (setq count (+ count 1)))))

        (and (= count 0)
	       ;; All pairs accounted for.
	       (point)))))))

; --- macro fortran-end-if-or-else ---
(defun fortran-end-if-or-else (&optional arg)
" Search forwards for first unmatched ENDIF or ELSE.  Return point or nil.
 Optional ARG=2 => if sitting at an ELSE, do not count it."

    (interactive "p")

  (let ((case-fold-search t))
    (if (save-excursion (beginning-of-line)
			(skip-chars-forward " \t0-9")
			(if (= arg 1)
			    (looking-at "end[ \t]*if\\b\\(else[ 	]*if\\b\\)\\|\\(else\\b\\)")
			  (looking-at "end[ \t]*if\\b")))
	;; Sitting on one.
	(match-beginning 0)
      ;; Search for one.  The point has been already been moved to first
      ;; letter on line but this should not cause troubles.
      (save-excursion
	(let ((count 1))
        (while (and (not (= count 0))
		      (not (eq (fortran-next-statement) 'last-statement))
		      ;; Keep local to subprogram.
		      (not (and (looking-at fortran-end-prog-re)
				(fortran-check-end-prog-re))))

	    (skip-chars-forward " \t0-9")
	    (cond

	     ((looking-at "end[ \t]*if\\b")
	      (setq count (- count 1)))

	     ((and (= count 1) (looking-at "\\(else[ 	]*if\\b\\)\\|\\(else\\b\\)"))
	      (setq count (- count 1)))

	     ((looking-at fortran-if-start-re)
	      (save-excursion
		(if (or
		     (looking-at ".*)[ \t]*then\\b[ \t]*[^ \t(=a-z0-9]")
		     (let (then-test) ; Multi-line if-then.
		       (while
			   (and (= (forward-line 1) 0)
				;; Search forward for then.
				(or (looking-at "     [^ 0\n]")
				    (looking-at "\t[1-9]"))
				(not
				 (setq then-test
				       (looking-at
					".*then\\b[ \t]*[^ \t(=a-z0-9]")))))
		       then-test))
		    (setq count (+ count 1)))))))

        (and (= count 0)
	     ;; All pairs accounted for.
	     (point)))))))

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

;(defun find-tagx (tagname &optional next other-window)
;  "Indentical to function FIND-TAG, except that find-tagx searches
;   the regular expression \\bTAGNAME\\b in the tags table,
;   which avoids accidental tags with longer names."

;  (interactive (if current-prefix-arg
;		   '(nil t)
;		 (find-tag-tag "Find tag: ")))
;  (let (buffer file linebeg startpos)
;    (save-excursion
;     (visit-tags-table-buffer)
;     (if (not next)
;	 (goto-char (point-min))
;       (setq tagname last-tag))
;     (setq last-tag tagname)
;     (while (progn
;	      (if (not (re-search-forward (concat "\\b" tagname "\\b") nil t))
;		  (error "No %sentries containing %s"
;			 (if next "more " "") tagname))
;	      (not (looking-at "[^\n\177]*\177"))))
;     (search-forward "\177")
;     (setq file (expand-file-name (file-of-tag)
;				  (file-name-directory tags-file-name)))
;     (setq linebeg
;	   (buffer-substring (1- (point))
;			     (save-excursion (beginning-of-line) (point))))
;     (search-forward ",")
;     (setq startpos (read (current-buffer))))
;    (if other-window
;	(find-file-other-window file)
;      (find-file file))
;    (widen)
;    (push-mark)
;    (let ((offset 1000)
;	  found
;	  (pat (concat "^" (regexp-quote linebeg))))
;      (or startpos (setq startpos (point-min)))
;      (while (and (not found)
;		  (progn
;		   (goto-char (- startpos offset))
;		   (not (bobp))))
;	(setq found
;	      (re-search-forward pat (+ startpos offset) t))
;	(setq offset (* 3 offset)))
;      (or found
;	  (re-search-forward pat nil t)
;	  (error "%s not found in %s" pat file)))
;    (beginning-of-line))
;  (setq tags-loop-form '(find-tag nil t))
;  ;; Return t in case used as the tags-loop-form.
;  t)

(defun is-tag (tagnam)
  " Returns t if tags table has tagnam, otherwise nil"

  (interactive "P")

  (save-excursion
    (visit-tags-table-buffer) (goto-char (point-min))
    (if (re-search-forward (concat "\\b" tagnam "\\b *(\177") nil t)
	(progn (beginning-of-line) t) nil)
))

; --- macro fortran-list-output-variables ---
(defun fortran-list-output-variables (&optional arg)
  " Lists variables set within a fortran subroutine"

   (interactive "P")

; --- Setup ---
   (let ((lastp (progn (end-of-fortran-subprogram) (point)))
	 (bnam (current-buffer)) strn pointnow)
(progn
     (pop-to-buffer "*vars*") (kill-buffer "*vars*")
     (pop-to-buffer "*vars*") (pop-to-buffer bnam)
     (fortran-first-executable-statement)
)

; --- Repeat while '=' within subroutine ---
     (while
	 (progn (re-search-forward "=" lastp 0) (< (point) lastp))

; ---  Don't bother if = part of commented lines or format statement
       (if (or (save-excursion
		 (beginning-of-line) (looking-at "^[Cc*]"))
	       (save-excursion
		 (fortran-previous-statement) (fortran-next-statement)
		 (looking-at "^......[ \t]*\\(format\\|write\\)")))
	   nil
	   (progn
	     (setq pointnow (point))
	     (backward-char 1) (fortran-backward-variable 1)
	     (setq strn
		   (buffer-substring
		     (point) (save-excursion (forward-sexp) (point))))

	     (set-buffer "*vars*") (goto-char (point-max))
	     (insert strn "\n")
	     (set-buffer bnam) (goto-char pointnow)
	     (message "%s" strn)
	     ))))
)

; --- macro fortran-forward-whitespace ---
(defun fortran-forward-whitespace (&optional arg)
" Forward past any f77 'white space', which includes blanks, tabs,
  comment lines and the first 6 columns in a continuation line.

 Optional ARG non-nil includes characters -+*/=, as whitespace, except
       if ARG=2 includes -+*/= (\",\" is excluded).
       if ARG=3 Comment and first 6 columns are not treated specially
                and moves past spaces, tabs, and newlines only.
                (just a 'normal' whitespace movement).
 Returns t if position point changes, otherwise nil."

   (interactive "P")
   (or (numberp arg) (setq arg (if arg 1 0)))


   (let ((pointnow (point)) wschar)
					;set arg=nil->0,  nonnil -> 1
     (setq wschar (cond ((= arg 0) " \t\n")
		         ((= arg 1) "-+*/=, \t\n")
		         ((= arg 2) "-+*/= \t\n")
		         ((= arg 3) (skip-chars-forward " \t\n") " \t\n")))

     (or (= arg 3)
     (while
	 (or (= (current-column) 5)
	     (looking-at (concat "[" wschar "]"))
	     (fortran-is-comment-line))
					; Skip over comments
       (if (fortran-is-comment-line) (fortran-next-statement))
					; and over whitespace
       (skip-chars-forward wschar)
					; and over the continuation char
       (if (= (current-column) 5) (forward-char 1))))
     (not (= (point) pointnow)))
)

; --- macro fortran-backward-whitespace ---
(defun fortran-backward-whitespace (&optional arg)
"  Backward past any f77 'white space', which includes blanks, tabs,
   comment lines and the first 6 columns in a continuation line.
   Optional ARG non-nil includes characters -+*/=, as whitespace.
   Returns t if position point changes, otherwise nil."

   (interactive "P")

   (let ((pointnow (point)) (wschar "-+*/=, \t\n"))
     (if (not arg) (setq wschar " \t\n"))
     (while
	 (or (= (current-column) 6)
	     (save-excursion (backward-char 1)
			     (looking-at (concat "[" wschar "]")))
	     (fortran-is-comment-line))
					; skip over comments
       (while (fortran-is-comment-line)
	 (previous-line 1) (end-of-line))
					; and over whitespace
       (skip-chars-backward wschar)
					; and over the continuation char
       (if (= (current-column) 6) (backward-char 1)))
     (not (= (point) pointnow)))
)

; --- macro fortran-forward-expr ---
(defun fortran-forward-expr (&optional arg arg2)
"Forward past one fortran expression.  Invokes fortran-forward-variable,
 and repeats for as long as a binop is encountered.
 Optional ARG repeats ARG times.
 Optional ARG2 controls flow over whitespace:
 Optional ARG2 controls what is meant by whitespace:
 ARG2 2 or 3  invokes fortran-forward-whitespace (arg2), else
 ARG2 nil     invokes fortran-forward-whitespace (1), else
 ARG2 non-nil invokes fortran-forward-whitespace (nil)."

   (interactive "p")
   (or arg (setq arg 1))
   (or (numberp arg2) (setq arg2 (if arg2 1 0)))

   (let ((relop "\\.\\(and\\|eq\\|g[et]\\|l[et]\\|ne\\|neqv\\|or\\)\\.")
	 (binop "-+*/") (nextwhite nil) (firstwhite arg2) psave)

     (if (= firstwhite 3) (setq nextwhite 3))

     (while (progn (setq arg (- arg 1)) (>= arg 0))
       (progn
	 ;; go forward one variable
	 (fortran-forward-variable 1 firstwhite)
	 ;; while looking-at binop|relop, jump past the next variable
	 (while (save-excursion
		  (fortran-forward-whitespace nextwhite)
		  (setq psave (point))
		  (cond ((looking-at (concat "[" binop "]")) (skip-chars-forward binop) (setq psave (point)) t)
			((looking-at relop) (skip-chars-forward relop) (setq psave (point)) t)))
	   (goto-char psave)
	   (fortran-forward-variable 1 (or nextwhite 1))))))
)

; --- macro fortran-forward-variable ---
(defun fortran-forward-variable (&optional arg arg2)
"Forward past one variable,array or function.
 Moves forward-whitespace; then moves past either:
 a word and subsequently an expression (..) if it follows word
 a string or concatenation of them
 Optional ARG repeats ARG times.

 Optional ARG2 controls what is meant by whitespace:
 ARG2 2 or 3  invokes fortran-forward-whitespace (arg2), else
 ARG2 nil     invokes fortran-forward-whitespace (1), else
 ARG2 non-nil invokes fortran-forward-whitespace (nil)."

   (interactive "p")
   (or arg (setq arg 1))
   (or (numberp arg2) (setq arg2 (if arg2 1 0)))


   (let ((vchar "A-Za-z0-9._%")
					; nexwhite: later fortran-forward-whitespace
	 (nextwhite nil) firstwhite)
					; arg to initial fortran-forward-whitespace
     (setq firstwhite (cond ((= arg2 0) 1)
			    ((= arg2 2) 2)
			    ((= arg2 3) (setq nextwhite 3))
			    (t nil)))
     (while (progn (setq arg (- arg 1)) (>= arg 0))
       (progn
	 (fortran-forward-whitespace firstwhite)
					; if at close of expr, go up one level
	 (if (looking-at ")") (progn (up-list 1)))
					; subsequent whitespace includes no binops
	 (let ((last (point)) (isquote nil) (unop  "-\\|\\.not\\.") psave)
	   ;; A fortran variable is either a string
	   (while (save-excursion
		    (fortran-forward-whitespace nextwhite)
		    (setq psave (point) isquote t)
		    (or (looking-at "'") (progn (setq isquote nil) (looking-at "//"))))
	     (progn (goto-char (1+ psave))
		    (skip-chars-forward (if isquote "^'" "^/"))
		    (forward-char 1)
					; recursively go forward to end of string
		    (or isquote (fortran-forward-variable 1 (or nextwhite 1)))
		    t))

	   ;; or it's a unop, name, array or number
	   (if (= (point) last)
	       (progn
					; a unop
		 (if (looking-at unop)
		     (progn (skip-chars-forward unop)
			    (fortran-forward-whitespace nextwhite)))
					; a name or number
		 (skip-chars-forward vchar)
					; jump over expression
		 (while (save-excursion
		       (fortran-forward-whitespace nextwhite)
		       (setq psave (point)) (looking-at "[(%]"))
		     (progn (goto-char psave) (forward-sexp) t))
		 ))))))
)

; --- macro fortran-backward-variable ---
(defun fortran-backward-variable (&optional arg)
"  Backward past one variable,array or function.
   Invokes fortran-backward-whitespace; then moves past a word
   and subsequently past a sexp if it precedes point.
   Optional ARG repeats ARG times."

   (interactive "p") (if (not arg) (setq arg 1))

   (let ((vchar "A-Za-z0-9._"))
     (while (> arg 0)
       (fortran-backward-whitespace 1)
       (if (save-excursion (backward-char 1) (looking-at "("))
	   (progn (up-list -1) (fortran-backward-whitespace 1)))
       (if (save-excursion (backward-char 1) (looking-at "[')]"))
	   (backward-sexp))
       (skip-chars-backward vchar)
       (setq arg (- arg 1))))
)

; --- macro fortran-next-call ---
(defun fortran-next-call (subend)
"   Find the next call to a fortran subprogram.  subend marks limit of search.
    Subexpression 1 (for e.g. match-beginning, match-end) is matched to function name."

;  (interactive)
;  (snit)

  (if (re-search-forward "^[1-9 \t].+[ \t]call +\\([a-zA-Z0-9_]+\\)" subend t)  ; Found a match
    (if (save-excursion (beginning-of-line) (not (looking-at ".*!.*call"))) ; If call is commented out look for another call
	t
        (fortran-next-call subend))))

; --- macro fortran-is-comment-line ---
(defun fortran-is-comment-line (&optional arg)
"  Returns t if point in a fortran comment, else nil"

   (interactive "p")
   (save-excursion (beginning-of-line) (looking-at "^[Cc*]"))
)

; --- macro compilation-fort-declare-fortran-variable-type ---
(defun compilation-fort-declare-fortran-variable-type ()
"  Get next undeclared variable in compilation buffer and insert
   a declaration in prior window"
   (interactive)

   (let (vnam version)
;  ... See if there are any missing declarations (SGI)
;    (snot)
;     (setq version
;     (if
;         (re-search-forward "Variable \"\\(\\w+\\)\" must be explicitly typed" (point-max) t)
;         6
;       (if
;  	 (re-search-forward "identifier \"\\(\\w+\\)\" has not been explicitly" (point-max) t)
;  	 7
; 	(if
; 	    (re-search-forward "explicit type must be specified for data object \"\\(\\w+\\)\"" (point-max) t)
; 	    8 nil))))

;  ... See if there are any missing declarations (dec linux fort)
   (setq version
   (if (re-search-forward "must have an explicit type[^^]+" (point-max) t) 6 nil))
   (if (= version 6)
       (progn (previous-line 1) (re-search-forward "\\([a-zA-Z0-9_]+\\)" (point-max) t)))

   (if (or (= version 6) (= version 7) (= version 8))
       (progn
	 (setq vnam (buffer-substring (match-beginning 1) (match-end 1)))
	 (setq bufnow  (buffer-name)
	       pointnow (point)
;	       buffer (find-file-noselect (car gud-last-last-frame))
	       window (display-buffer bufnow))

	 (other-window -1) (fortran-first-executable-statement)
	 (if (or (and (>= (string-to-char vnam) ?I) (<= (string-to-char vnam) ?O))
		 (and (>= (string-to-char vnam) ?i) (<= (string-to-char vnam) ?o)))
	     (progn (re-search-backward "^ +integer")
		    (if (looking-at "^ +integer w(") (re-search-backward "^ +integer")))
	   (re-search-backward "^ +double precision"))
	 (fortran-next-statement)
	 (while (progn (previous-line 1) (looking-at "^[Cc]\\|^$")))
	 (end-of-line)
	 (insert-string (concat "," vnam))
	 (backward-word 1) (downcase-word 1)
	 (other-window 1)
	 (message "inserting variable %s" vnam))
	 (message "no variables found"))
))

; --- macro f77-to-quad ---
(defun f77-to-quad ()
  "Converts intrinsic functions to quad for compilation with
   arbitrary precision."
(interactive)
   (save-excursion
;    (untabify (point-min) (point-max))

     (beginning-of-buffer)
     (while (progn (while (looking-at "^[Cc*]") (next-line 1))
		   (< (point) (point-max)))
; skip over comments
       (while (looking-at "^[Cc*]") (next-line 1))
; skip over data and format statements
       (while (looking-at "^...... *\\(format\\|data\\)") (fortran-next-statement))
       (fortran-const-to-quad (point) (progn (end-of-line) (point)) nil)
       (forward-char 1)
     )


     (beginning-of-buffer)
     (while (progn (while (looking-at "^[Cc*]") (next-line 1))
		   (< (point) (point-max)))
; skip over comments
       (while (looking-at "^[Cc*]") (next-line 1))
; skip over format statements
       (while (looking-at "^...... *format") (fortran-next-statement))
       (fortran-intrsc-to-quad (point) (progn (end-of-line) (point)) nil)
       (forward-char 1)
     )
))


; --- macro fortran-cnvt-buffer-to-quad ---
(defun fortran-cnvt-buffer-to-quad ()
  "Converts intrinsic functions to quad for compilation with
   arbitrary precision."
(interactive)
   (save-excursion
     (untabify (point-min) (point-max))
     (beginning-of-buffer)
     (while (< (point) (point-max))
; skip over comments
        (while (looking-at "^[Cc*]") (next-line 1))
        (fortran-intrsc-to-quad (point) (progn (end-of-line) (point)) t)
        (forward-char 1)
     )
))

(defun fortran-const-to-quad (startposn posn que)
  "converts fortran constants to quad precision in region (STARTPOSN, POSN).
   If QUE (third arg) is t, asks interactively.
   Returns change in size of region."
  (interactive)
					; Set up for query replace
  (if que (setq query ??) (setq query ?a))
  (let ((pend posn))
					; replace isolated numbers beginning with a digit or "."
	(goto-char startposn)
	(while (and (< (point) pend)
		    (re-search-forward "\\b\\([0-9][0-9 ]*\\.*[0-9 ]*d *[+-]?[0-9]+\\)\\|\\b\\(\\.[0-9][0-9 ]*d *[+-]?[0-9]+\\)" pend t))
	  (progn
	    (setq query (get-query "Convert?" query))
	    (if (or (char-equal query ?a) (char-equal query ?y))
		(progn (replace-match "dble(\\1\\2)" t)
		       (setq pend (+ pend 6))))))


				; replace numbers following .log., beginning with a digit, eg == 2d0
    (goto-char startposn)
    (while (and (< (point) pend)
		(re-search-forward "\\.\\(and\\|eq\\|false\\|ne\\|not\\|or\\|true\\)\\.\\([0-9][0-9 ]*\\.*[0-9 ]*d *[+-]?[0-9]+\\)" pend t))
      (progn
	(setq query (get-query "Convert?" query))
	(if (or (char-equal query ?a) (char-equal query ?y))
	    (progn (replace-match ".\\1.dble(\\2)" t)
		   (setq pend (+ pend 6))))))

					; handle cases eg == .2d0
;      (goto-char startposn)
;      (while (re-search-forward "\\.\\(and\\|eq\\|false\\|ne\\|not\\|or\\|true\\)\\.\\(\\.*[0-9][0-9 ]*d *[+-]?[0-9]+\\)" pend t)
;	(progn
;	  (setq query (get-query "Convert?" query))
;	  (if (or (char-equal query ?a) (char-equal query ?y))
;	      (progn (replace-match ".\\1.dble(\\2)" t)
;		     (setq pend (+ pend 6))))))

    (goto-char pend)
    (- pend posn)))

(defun fortran-intrsc-to-quad (startposn posn que)
  "converts intrinsic fortran functions to quad in region (STARTPOSN, POSN).
   If QUE (third arg) is t, asks interactively.
   Returns change in size of region."
  (interactive)
; Set up for query replace
   (if que (setq query ??) (setq query ?a))
   (setq savpos posn)

; Handle all cases where dfnct -> qfnct
   (goto-char startposn)
   (while (re-search-forward
           "\\bd\\(imag\\|cmplx\\|conjg\\|log\\|log10\\)\\b" posn t) (progn
      (setq query (get-query "Convert?" query))
      (if (or (char-equal query ?a) (char-equal query ?y))
         (progn (backward-word 1) (delete-char 1) (insert-string "q"))
   )))


; Handle all cases where cdfnct -> fnct
    (goto-char startposn)
    (while (re-search-forward
            "\\bcd\\(abs\\|cos\\|exp\\|sin\\)\\b"  posn t) (progn
      (setq query (get-query "Convert?" query))
      (if (or (char-equal query ?a) (char-equal query ?y))
         (progn (backward-word 1) (delete-char 2)
              (setq posn (- posn 2))
   ))))


; Handle all cases where cdfnct -> cqfnct
    (goto-char startposn)
    (while (re-search-forward
            "\\bcd\\(log\\)\\b"  posn t) (progn
      (setq query (get-query "Convert?" query))
      (if (or (char-equal query ?a) (char-equal query ?y))
         (progn (backward-word 1) (delete-char 2) (insert-string "cq")
              (setq posn (- posn 0))
   ))))


; Handle dmin
    (goto-char startposn)
    (while (re-search-forward
            "\\b\\(dmin1\\|dmax1\\)\\b"  posn t) (progn
      (setq query (get-query "Convert?" query))
      (if (or (char-equal query ?a) (char-equal query ?y))
	  (progn (backward-char 1) (delete-char 1) (backward-word 1)
		  (delete-char 1) (insert-string "x") (setq posn (- posn 2))
   ))))

; Handle nint,min,max
;    (goto-char startposn)
;    (while (re-search-forward
;            "\\b\\(nint\\|min\\|max\\)\\b"  posn t) (progn
;      (setq query (get-query "Convert?" query))
;      (if (or (char-equal query ?a) (char-equal query ?y))
;         (progn (backward-word 1) (insert-string "x")
;              (setq posn (+ posn 1))
;   ))))

; Handle all cases where dfnct -> fnct
   (goto-char startposn)
   (while (re-search-forward
           "\\bd\\(abs\\|acos\\|asin\\|atan\\|atan2\\|cos\\|cosh\\|exp\\|sign\\|sin\\|sinh\\|sqrt\\|tan\\|tanh\\)\\b"  posn t) (progn
      (setq query (get-query "Convert?" query))
      (if (or (char-equal query ?a) (char-equal query ?y))
         (progn (backward-word 1) (delete-char 1)
            (setq posn (- posn 1))
   ))))

 (- posn savpos))                    ; Return with change in region size

(defun fortran-replace-do-enddo()
" Prompts for a logical unit, then replaces do-enddo with f77 standard.
   Fortran-do-enddo block is marked by fortran-mark-do."
  (interactive)
  (fortran-mark-do)
  (let (numstr                         ;line numstr to use
	(default "30"))                 ;default str -- not used now

	(setq numstr (read-from-minibuffer
		      (format "Use line number for  '%s' :"
			      (buffer-substring (point) (save-excursion (end-of-line) (point))))
		      nil nil nil nil))
	(forward-char 2)
	(delete-horizontal-space)
	(insert-string (concat "  " numstr "  "))
	(exchange-point-and-mark)
	(beginning-of-line)
	(insert-string (concat numstr "     continue"))
	(kill-line nil)
	(fortran-indent-line)))

(defun fortran-pad-buffer()
" Pads binops, do-loops and relops spaces in buffer."
  (interactive)
  (goto-char (point-min))  (fortran-pad-relop)
  (goto-char (point-min))  (fortran-pad-binop)
  (goto-char (point-min))  (fortran-pad-do-loops)
)

(defun fortran-pad-relop()
" Pads .op. with spaces using query-replace-regexp."
  (interactive)
  (let ((comment-chars "c!*")
	(fortran-relops
;       ("and" "or" "not" "lt" "le" "eq" "ge" "gt" "ne")
         "and\\|eq\\|g[et]\\|l[et]\\|n\\(e\\|ot\\)\\|or"))

    (query-replace-regexp
     (concat "\\(^ .+[^ ]\\) *\\.\\(" fortran-relops "\\)\\. *")
     "\\1 .\\2. " nil)))

(defun fortran-pad-binop()
" Pads =,+,- with spaces using query-replace-regexp."
  (interactive)
  (save-excursion
    (query-replace-regexp "\\(^ .+[^ ]\\) *\\([=]\\) *" "\\1 \\2 " nil))

  (save-excursion
  (query-replace-regexp "\\(^ .+[^ ]\\) *\\([+]\\) *" "\\1 \\2 " nil))

  (save-excursion
  (query-replace-regexp "\\(^ .+[^ ]\\) *\\([-]\\) *" "\\1 \\2 " nil)))


(defun fortran-pad-do-loops()
" Pads do-loops with spaces using query-replace-regexp."
  (interactive)
  (query-replace-regexp "do +\\(\\w+\\) +\\(\\w+\\) *= *\\(\\w+\\) *, *\\(\\w+\\)" "do  \\1  \\2 = \\3, \\4" nil)
  (query-replace-regexp "do +\\(\\w+\\) *= *\\(\\w+\\) *, *\\(\\w+\\)" "do  \\1 = \\2, \\3" nil)
)

(defun fortran-ftnchek-clean-unused-var()
" Uses output from ftnchek to eliminate unused variables.  Run ftnchek in
 a shell and put the cursor before 'File'.  Uses remember-point"
  (interactive)

  (search-forward "Variables declared but never referenced")
  (search-backward "\"")
  (beginning-of-line)
  (search-forward "\"")
  (let (vnam
	(fname (buffer-substring
		(point) (progn (search-forward "\"") (- (point) 1))))
	(module (buffer-substring
		 (progn (search-forward "module ") (point))
		 (progn (search-forward ":") (- (point) 1))))
	(pmax (save-excursion (re-search-forward "^$") (point))))
    (search-forward "Variables declared but never referenced")
    (beginning-of-line) (forward-line 1)
    (find-file-other-window fname) (goto-char (point-min))
    (re-search-forward (concat "^ .+\\(subroutine\\|function\\) +" module "\\b"))
    (remember-point)
    (other-window -1)
    (while (and (re-search-forward "[A-Za-z]") (< (point) pmax))
      (backward-char 1)
      (setq vnam (buffer-substring
		(point) (progn (re-search-forward "[* \t]") (- (point) 1))))
      (message "name is %s" vnam)
      (other-window 1) (recover-point)
      (query-replace-regexp (concat "\\b" vnam "\\b,*") "")
      (other-window -1))
)
)

(defun fortran-fort-clean-unused-var()
" Uses output from DEC fort to eliminate unused variables.  Run fort in
 a shell and put the cursor before 'File'."
  (interactive)

; (snot)

  (let (vnam lineno
        (fname (buffer-substring
                (progn (re-search-forward "info: ") (point))
                (progn (re-search-forward ",") (- (point) 1))))
)

    (find-file-other-window fname)
;    (re-search-forward (concat "^ .+\\(subroutine\\|function\\) +" module "\\b"))
;    (remember-point)
    (other-window -1)

    (while (re-search-forward "This variable has not been used" (point-max) t)
      (re-search-forward "\\[\\([^]]+\\)")
      (setq vnam (downcase (buffer-substring (match-beginning 1) (match-end 1))))
      (message "name is %s" vnam)
      (beginning-of-line)
      (setq lineno (progn (re-search-forward "line +\\([0-9]+\\)") (match-string 1)))
      (forward-line 1)
      (other-window 1) (goto-line lineno) (beginning-of-line)
      (query-replace-regexp (concat ", *\\b" vnam "\\b\\(([^)]+)\\)*\\|\\b" vnam "\\b\\(([^)]+)\\)* *,\\|\\b" vnam "\\b\\(([^)]+)\\)*") "")
      (other-window -1))
))


;; (defun fortran-flint-clean-unused-var()
;; " Uses output from SGI ftnlint to eliminate unused variables.  Run flint in
;;  a shell and put the cursor before 'File'.  Uses remember-point"
;;   (interactive)

;; ; (snot)

;;   (let (vnam
;;         (fname (buffer-substring
;;                 (progn (re-search-forward "ftnlint  +") (point))
;;                 (progn (re-search-forward "$") (- (point) 0))))
;;         (module (buffer-substring
;;                   (progn (search-forward "ftnlint message") (beginning-of-line) (skip-chars-forward " \t") (point))
;;                    (progn (search-forward "(") (- (point) 1))))
;;         (pmax (save-excursion (re-search-forward "ftnlint  ") (beginning-of-line) (point))))

;;     (find-file-other-window fname) (goto-char (point-min))
;;     (re-search-forward (concat "^ .+\\(subroutine\\|function\\) +" module "\\b"))
;;     (remember-point)
;;     (other-window -1)

;;     (while (re-search-forward "never used\\|No references to parameter" pmax t)
;;       (beginning-of-line) (search-forward "\"")
;;       (setq vnam (buffer-substring
;;                    (point) (progn (search-forward "\"") (- (point) 1))))
;;       (message "name is %s" vnam)
;;       (forward-line 1)
;;       (other-window 1) (recover-point)
;;       (query-replace-regexp (concat ", *\\b" vnam "\\b\\(([^)]+)\\)*\\|\\b" vnam "\\b\\(([^)]+)\\)* *,\\|\\b" vnam "\\b\\(([^)]+)\\)*") "")
;;       (other-window -1))
;; ))

;; (defun fortran-flint-find-undef-var()
;; " Uses output from flint to eliminate unused variables.  Run flint in
;;  a shell and put the cursor before 'File'.  Uses remember-point"
;;   (interactive)

;;   (let (vnam
;;         (fname (buffer-substring
;;                 (progn (re-search-forward "ftnlint  +") (point))
;;                 (progn (re-search-forward "$") (- (point) 0))))
;;         (module (buffer-substring
;;                   (progn (search-forward "ftnlint message") (beginning-of-line) (skip-chars-forward " \t") (point))
;;                    (progn (search-forward "(") (- (point) 1))))
;;         (pmax (save-excursion (re-search-forward "ftnlint  ") (beginning-of-line) (point))))

;;     (find-file-other-window fname) (goto-char (point-min))
;;     (re-search-forward (concat "^ .+\\(subroutine\\|function\\) +" module "\\b"))
;;     (remember-point)
;;     (other-window -1)

;; ;(snot)

;;     (while (re-search-forward "never assigned a value\\|may be used before it is assigned" pmax t)
;;       (beginning-of-line) (search-forward "\"")
;;       (setq vnam (buffer-substring
;;                    (point) (progn (search-forward "\"") (- (point) 1))))
;;       (message "name is %s" vnam)
;;       (forward-line 1)
;;       (other-window 1) (recover-point) (fortran-first-executable-statement)
;;       (query-replace-regexp (concat ", *\\b" vnam "\\b\\(([^)]+)\\)*\\|\\b" vnam "\\b\\(([^)]+)\\)* *,\\|\\b" vnam "\\b\\(([^)]+)\\)*") "")
;;       (other-window -1))
;; ))


(defun fortran-extract-struc-list (&optional arg)
"Get structure list."

   (interactive "p")
   (let ((bnam (current-buffer)) vlist)

     (re-search-backward "^Cr +off +offe")
     (next-line 1)
     (setq vlist
	   (buffer-substring
	    (point) (progn (fortran-next-statement) (point))))
     (set-buffer (get-buffer-create " *snot*"))
     (erase-buffer)
     (insert-string vlist)
     (goto-char 1)
     (while (< (point) (point-max))
       (delete-char 4)
       (if (looking-at " ") (delete-char 1))
       (if (looking-at " ") (delete-char 1))
       (if (looking-at " ") (delete-char 1))
       (if (looking-at "[0-9]")
	   (progn (skip-chars-forward "0-9") (insert-string ",")))
       (kill-line)
       (delete-char 1)
       )

     (setq vlist
	   (buffer-substring
	    (point-min) (point-max)))

     (kill-buffer (current-buffer))
     (set-buffer bnam)
     (insert-string (concat "      data ilists /" vlist "/\n"))
     (backward-char 2))
)

; --- macro fortran-sort-declaration-list ---
(defun fortran-sort-declaration-list (&optional arg)

   (interactive)

   (let ((fortran-type-types
	  (concat "byte\\|c\\(haracter\\|om\\(mon\\|plex\\)\\)\\|"
                  "d\\(ata\\|imension\\|ouble"
                  "[ \t]*\\(complex\\|precision\\)\\)\\|"
                  "e\\(nd[ \t]*\\(map\\|structure\\|union\\)\\|"
                  "quivalence\\|xternal\\)\\|"
                  "i\\(mplicit[ \t]*\\(byte\\|"
                  "c\\(haracter\\|omplex\\)\\|"
                  "double[ \t]*\\(complex\\|precision\\)\\|"
                  "integer\\|logical\\|none\\|real\\)\\|"
                  "nt\\(eger\\|rinsic\\)\\)\\|"
                  "logical\\|map\\|none\\|parameter\\|re\\(al\\|cord\\)\\|"
                  "s\\(ave\\|tructure\\)\\|union"))
	 subrstrng eod
	 (thisbuf (current-buffer)))

; .. Find start of line
     (beginning-of-line)
     (fortran-next-statement)
     (fortran-backward-whitespace)
     (setq eod (point))
     (fortran-next-statement)
     (fortran-previous-statement)
; .. Function body, conditional on being at a declaration
     (if (looking-at (concat "^[ \t] +" "\\<\\(" fortran-type-types "\\)\\>"))
	 (progn
	   (re-search-forward (concat "^[ \t] +" "\\<\\(" fortran-type-types "\\)\\>"))
	   (fortran-forward-whitespace)
	   (setq subrstrng (buffer-substring (point) eod))
	   (pop-to-buffer "*compare*") (erase-buffer) (insert "     ." subrstrng "\n")
	   (goto-char (point-min))
					; put each declaration into a separate column
	   (while (< (point) (point-max))
	     (fortran-forward-variable 1) (insert ",") (fortran-forward-whitespace) (insert "\n     .")
	     (if (looking-at "[ ,]") (delete-char 1))
	     (while (and (< (point) (point-max)) (looking-at "$")) (delete-char 7) (while (looking-at " ") (delete-char 1))))
	   (sort-lines nil (point-min) (point-max))
	   (goto-char (point-min))
	   (while (looking-at "^$") (delete-char 1))
	   (while (looking-at "^     \\.$") (delete-char 7))
	   (pop-to-buffer thisbuf)
	   (setq subrstrng (buffer-substring (save-excursion (beginning-of-line) (point)) (point)))
	   (fortran-next-statement)
	   (fortran-previous-statement)
	   (insert subrstrng "\n")
	   (pop-to-buffer "*compare*")
	   (setq subrstrng (buffer-substring (point-min) (point-max)))
	   (pop-to-buffer thisbuf)
	   (insert subrstrng "\n")
	   )
       (message "Not at a declaration --- cannot sort"))
))

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

(load "add-log")
(defun add-doc-iso8601-time-string ()
  (if change-log-time-zone-rule
      (let ((tz (getenv "TZ"))
	    (now (current-time)))
	(unwind-protect
	    (progn
	      (set-time-zone-rule
	       change-log-time-zone-rule)
	      (concat
	       (format-time-string "%Y-%m-%d " now)
	       (add-log-iso8601-time-zone now)))
	  (set-time-zone-rule tz)))
    (format-time-string "%d %h %y")))

; --- macro fortran-split-dimension-statement ---
(defun fortran-split-dimension-statement ()
"   Splits fortran 'dimension' statement into separate integer and
    double precision declarations"
   (interactive)

   (let ((strnint "") (strndbl ""))

     (fortran-previous-statement)
     (fortran-next-statement)
     (beginning-of-line)
     (insert (buffer-substring
	      (point) (save-excursion (fortran-next-statement) (point))))
     (fortran-comment-region
      (point) (save-excursion (fortran-next-statement) (point)) nil)
     (fortran-previous-statement)

     (while (progn (end-of-line) (looking-at "\n     [^ ]"))
       (progn (delete-char 7) (delete-horizontal-space)))

; --- break dimension declaration into integer and double precision
     (while (progn (beginning-of-line) (looking-at "^ +dimension +[a-zA-Z]"))
       (if (looking-at "^ +dimension [A-HP-Za-hp-z]")
	   (progn (re-search-forward "^ +dimension +")
		  (setq strndbl (concat strndbl
					(buffer-substring (point) (save-excursion (fortran-forward-variable 1) (point)))
					","))
		  (fortran-kill-word 1)
		  (delete-horizontal-space) (insert " ")
		  (if (looking-at ",") (delete-char 1))
		  )
	   (progn (re-search-forward "^ +dimension +")
		  (setq strnint (concat strnint
					(buffer-substring (point) (save-excursion (fortran-forward-variable 1) (point)))
					","))
		  (fortran-kill-word 1)
		  (delete-horizontal-space) (insert " ")
		  (if (looking-at ",") (delete-char 1))
		  )
	 ))

; --- delete now-empty dimension statement
     (beginning-of-line) (kill-line 1) (fortran-next-statement)

;	  (snit)

; --- create integer and double precision declarations
   (if (string= strnint "") nil
     (progn
;       (fortran-next-statement)
       (insert (concat "      integer " strnint))
       (delete-backward-char 1)
       (insert "\n")))
   (if (string= strndbl "") nil
     (progn
;       (fortran-next-statement)
       (insert (concat "      double precision " strndbl))
       (delete-backward-char 1)
       (insert "\n")))

; --- delete blank lines following current line
     (while (progn (beginning-of-line) (looking-at " *$"))
       (progn (kill-line 1)))

))

; --- macro fortran-msm-to-standard ---
(defun fortran-msm-to-standard ()
"   Partially convert subroutine from MSM to standard style"
   (interactive)

 (let (
       (eor (progn (end-of-fortran-subprogram) (point)))
       sblank scomment eos int)

   (fortran-pad-buffer)

   ; Separate 'dimension' statement into 'integer' and 'dp'
   (beginning-of-fortran-subprogram) (fortran-next-statement)
   (if (re-search-forward "^ +dimension" eor t)
       (progn
         (replace-match "      double precision")

         ; create integer declaration ; int points to end of it
         (save-excursion
           (beginning-of-line)
           (insert-string "      integer \n") (setq int (- (point) 1)))

         ; move integer arrays to integer declaration
         (while
             (progn
               ; find eos=end-of-statement
               (save-excursion
                 (beginning-of-line) (fortran-next-statement)  (setq eos (- (point) 1)))
               ; do until past end-of-statement
               (fortran-forward-whitespace)
               (< (point) eos))

           (if (looking-at "[A-Ha-hP-Zp-z]")
               (progn (fortran-forward-variable 1) (if (looking-at ",") (forward-char 1)))
             (fortran-kill-word 1)
             (save-excursion
               (goto-char int) (yank) (insert-string ",") (setq int (point)))
             (if (looking-at ",") (delete-char 1))))

         ; cleanup
         (progn
           (goto-char int)
           (backward-char 1)
           (if (looking-at ",") (delete-char 1)))

    ))

   ; Adjust comment lines
   (fortran-first-executable-statement)
   (while (re-search-forward "^c *\\(---\\|\\.\\.\\.\\)\\( *\\) " eor t)
     (progn
       (setq scomment (buffer-substring (match-beginning 1) (match-end 1)))
       (setq sblank (buffer-substring (match-beginning 2) (match-end 2)))
       (if (string= scomment "---")
          (progn (beginning-of-line) (capitalize-word 2) (end-of-line) (insert-string " ---"))
        (progn (beginning-of-line) (capitalize-word 3)))
       (replace-match "C \\2\\1 ")))

   (beginning-of-fortran-subprogram)
   (query-replace-regexp "s_atom" "ssite")
   (beginning-of-fortran-subprogram)
   (query-replace-regexp "s_" "s")

   (beginning-of-fortran-subprogram) (fortran-next-statement)
   (insert "      implicit none\nC ... Passed parameters\nC ... Local parameters\n")
   (insert-string "      integer stdo,nglob\n")

   (fortran-first-executable-statement)
   (query-replace-regexp "write *(6," "write(stdo,")


   (fortran-first-executable-statement)
   (insert "C ... Heap\n      integer w(1)\n      common /w/ w\n")
   (insert-string "      stdo = nglob('stdo')\n")

   (fortran-first-executable-statement)
   (replace-regexp "dreal *(" "dble(")
   (beginning-of-fortran-subprogram)
   (replace-regexp "complex\\*16" "double complex")
))

; --- macro fortran-replace-do-with-do-enddo ---
(defun fortran-replace-do-with-do-enddo()
" Replaces f77 standard fortran do with do-enddo"
  (interactive)
  (beginning-of-line)
; Point at start of do loop
  (if (fortran-check-for-matching-do) (fortran-previous-block))
; Only perform replacement if start and end can be found
  (if (looking-at "[ \t]+do[ \t]*\\([0-9]+\\)") (progn
     (fortran-next-block)
     (if (fortran-check-for-matching-do) (progn
;    Perform replacement with safe assumption startloop and endloop will be correct
     (let (
        strn
        pointnow pointsave
        (endloop (point))
	(startloop (progn (fortran-previous-block) (point))))

					;patch enddo
	(goto-char endloop)
	(if (not (looking-at "[ \t]+[0-9]+[ \t]+continue"))
		 (progn
		   (setq pointnow (point))
		   (setq strn (buffer-substring
			    (point) (save-excursion (skip-chars-forward "[ \t]+[0-9]+[ \t]+") (point))))
		   (fortran-next-statement)
		   (beginning-of-line) (setq endloop (point))
		   (insert (concat strn "continue\n"))
		   (goto-char pointnow)
		   (skip-chars-forward "[ \t]+")
                   (while (looking-at "[0-9]") (delete-char 1) (insert " "))))

	(goto-char endloop)
	(beginning-of-line) (fortran-next-statement) (insert "      enddo\n")

					;patch "do"
	(goto-char startloop)
	(re-search-forward "\\([ \t]+do[ \t]*\\)\\([0-9]+\\)")
	(replace-match "\\1" t nil nil nil)
	(delete-horizontal-space)
	(insert-string "  ")

;    End of replacement
     )
)))))

; --- macro partks-to-gtv ---
(defun partks-to-gtv ()
" Approximately convert partks call to gtv call."
    (interactive)

  (let ((pointnow (point))
	(thisbuf (current-buffer))
	(ctrlbuf "m_rdctrl.f")
	catnam cast toknam varnam helpstr defval strucnam)

    (save-excursion
      (re-search-backward "getcat")
      (goto-char (match-end 0)) (down-list 1) (fortran-forward-variable 1)
      (re-search-forward "'")
      (setq catnam (buffer-substring (point) (progn (skip-chars-forward "[A-Za-z0-9]+") (point)))))

    (fortran-next-statement) (fortran-previous-statement)
    (re-search-forward "partks") (down-list 1) (fortran-forward-variable 1)


    (save-excursion
      (re-search-forward "'")
      (setq toknam (buffer-substring (point) (progn (skip-chars-forward "[A-Za-z0-9]+") (point)))))
    (fortran-forward-variable 2)
    (re-search-forward "'")
    (setq strucnam (buffer-substring (point) (progn (skip-chars-forward "[A-Za-z0-9]+") (point))))
    (re-search-forward " +")
    (setq varnam (buffer-substring (point) (progn (skip-chars-forward "[A-Za-z0-9_,]+") (point))))
    (re-search-forward ",")
					;   For SITE and SPEC
    (if (or (looking-at "ssite") (looking-at "sspec"))
	    (setq catnam (concat catnam "_ATOM")))
    (fortran-forward-variable 2) (re-search-forward "'") (backward-char 1)
    (setq helpstr (buffer-substring (point) (progn (fortran-forward-variable 1) (point))))
    (fortran-forward-variable 2) (fortran-backward-variable 1)
    (if (looking-at "4") (setq cast 4))
    (if (looking-at "2") (setq cast 2))
    (if (looking-at "1") (setq cast 1))
    (if (looking-at "0") (setq cast 0))

    (if (= cast 4) (setq defval "def_r8=0d0"))
    (if (= cast 2) (setq defval "def_i4=0"))
    (if (= cast 1) (setq defval "nmin=10"))
    (if (= cast 0) (setq defval "def_lg=F"))

    (pop-to-buffer ctrlbuf)

    (insert-string (concat
		    "      "
		    "nm='" catnam "_" toknam "'; call gtv(trim(nm),tksw(prgn,nm)," strucnam "_" varnam ",\n"
		    "     .  " defval ",note=" helpstr ")\n"))
))

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



; (cancel-debug-on-entry 'fortran-msm-to-standard)
; (local-set-key "\C-cm" 'fortran-msm-to-standard)

; (defun snot () (interactive "P") )
; (debug-on-entry 'snot)


(provide 'fortran)

;;; fortran.el ends here


