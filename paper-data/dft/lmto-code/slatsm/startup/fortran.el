;;; fortran.el --- Fortran mode for GNU Emacs

;;; Copyright (c) 1986, 1993, 1994, 1995 Free Software Foundation, Inc.

;; Author: Michael D. Prange <prange@erl.mit.edu>
;; Maintainer: bug-fortran-mode@erl.mit.edu
;; Version 1.30.6 (July 27, 1995)
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
;; along with GNU Emacs; see the file COPYING.  If not, write to
;; the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.

;;; Commentary:

;; Fortran mode has been upgraded and is now maintained by Stephen A. Wood
;; (saw@cebaf.gov).  It now will use either fixed format continuation line
;; markers (character in 6th column), or tab format continuation line style
;; (digit after a TAB character.)  A auto-fill mode has been added to
;; automatically wrap fortran lines that get too long.

;; We acknowledge many contributions and valuable suggestions by
;; Lawrence R. Dodd, Ralf Fassel, Ralph Finch, Stephen Gildea,
;; Dr. Anil Gokhale, Ulrich Mueller, Mark Neale, Eric Prestemon, 
;; Gary Sabot and Richard Stallman.

;;; This file may be used with GNU Emacs version 18.xx if the following
;;; variable and function substitutions are made.
;;;  Replace:
;;;   frame-width                           with screen-width
;;;   auto-fill-function                    with auto-fill-hook
;;;   comment-indent-function               with comment-indent-hook
;;;   (setq unread-command-events (list c)) with (setq unread-command-char c)

;;; Bugs to bug-fortran-mode@erl.mit.edu

(defconst fortran-mode-version "version 1.30.6")

;;; Code:

;;;###autoload
(defvar fortran-tab-mode-default nil
  "*Default tabbing/carriage control style for empty files in Fortran mode.
A value of t specifies tab-digit style of continuation control.
A value of nil specifies that continuation lines are marked
with a character in column 6.")

;; Buffer local, used to display mode line.
(defvar fortran-tab-mode-string nil
  "String to appear in mode line when TAB format mode is on.")

(defvar fortran-do-indent 3
  "*Extra indentation applied to DO blocks.")

(defvar fortran-if-indent 3
  "*Extra indentation applied to IF blocks.")

(defvar fortran-structure-indent 3
  "*Extra indentation applied to STRUCTURE, UNION, MAP and INTERFACE blocks.")

(defvar fortran-continuation-indent 5
  "*Extra indentation applied to Fortran continuation lines.")

(defvar fortran-comment-indent-style 'fixed
  "*nil forces comment lines not to be touched,
'fixed makes fixed comment indentation to `fortran-comment-line-extra-indent'
columns beyond `fortran-minimum-statement-indent-fixed' (for
`indent-tabs-mode' of nil) or `fortran-minimum-statement-indent-tab' (for
`indent-tabs-mode' of t), and 'relative indents to current
Fortran indentation plus `fortran-comment-line-extra-indent'.")

(defvar fortran-comment-line-extra-indent 0
  "*Amount of extra indentation for text within full-line comments.")

(defvar comment-line-start nil
  "*Delimiter inserted to start new full-line comment.")

(defvar comment-line-start-skip nil
  "*Regexp to match the start of a full-line comment.")

(defvar fortran-minimum-statement-indent-fixed 6
  "*Minimum statement indentation for fixed format continuation style.")

(defvar fortran-minimum-statement-indent-tab (max tab-width 6)
  "*Minimum statement indentation for TAB format continuation style.")

;; Note that this is documented in the v18 manuals as being a string
;; of length one rather than a single character.
;; The code in this file accepts either format for compatibility.
(defvar fortran-comment-indent-char " "
  "*Single-character string inserted for Fortran comment indentation.
Normally a space.")

(defvar fortran-line-number-indent 1
  "*Maximum indentation for Fortran line numbers.
5 means right-justify them within their five-column field.")

(defvar fortran-check-all-num-for-matching-do nil
  "*Non-nil causes all numbered lines to be treated as possible DO loop ends.")

(defvar fortran-blink-matching-if nil
  "*Non-nil causes \\[fortran-indent-line] on ENDIF statement to blink on matching IF.
Also, from an ENDDO statement blink on matching DO [WHILE] statement.")

(defvar fortran-continuation-string "$"
  "*Single-character string used for Fortran continuation lines.
In fixed format continuation style, this character is inserted in
column 6 by \\[fortran-split-line] to begin a continuation line.
Also, if \\[fortran-indent-line] finds this at the beginning of a line, it will
convert the line into a continuation line of the appropriate style.
Normally $.")

(defvar fortran-comment-region "c$$$"
  "*String inserted by \\[fortran-comment-region]\
 at start of each line in region.")

(defvar fortran-electric-line-number t
  "*Non-nil causes line number digits to be moved to the correct column as\
 typed.")

(defvar fortran-startup-message t
  "*Non-nil displays a startup message when Fortran mode is first called.")

(defvar fortran-column-ruler-fixed
  "0   4 6  10        20        30        40        5\
\0        60        70\n\
\[   ]|{   |    |    |    |    |    |    |    |    \
\|    |    |    |    |}\n"
  "*String displayed above current line by \\[fortran-column-ruler].
This variable used in fixed format mode.")

(defvar fortran-column-ruler-tab
  "0       810        20        30        40        5\
\0        60        70\n\
\[   ]|  { |    |    |    |    |    |    |    |    \
\|    |    |    |    |}\n"
  "*String displayed above current line by \\[fortran-column-ruler].
This variable used in TAB format mode.")

(defconst bug-fortran-mode "bug-fortran-mode@erl.mit.edu"
  "Address of mailing list for Fortran mode bugs.")

(defvar fortran-mode-syntax-table nil
  "Syntax table in use in Fortran mode buffers.")

(defvar fortran-analyze-depth 100
  "Number of lines to scan to determine whether to use fixed or TAB format\
 style.")

(defvar fortran-break-before-delimiters t
  "*Non-nil causes `fortran-fill' to break lines before delimiters.")

(if fortran-mode-syntax-table
    ()
  (setq fortran-mode-syntax-table (make-syntax-table))
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
  (modify-syntax-entry ?. "w" fortran-mode-syntax-table)
  (modify-syntax-entry ?_ "w" fortran-mode-syntax-table)
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

(let ((comment-chars "c!*")
      (fortran-type-types
;     (make-regexp
;      (let ((simple-types '("character" "byte" "integer" "logical"
;			    "none" "real" "complex"
;			    "double[ \t]*precision" "double[ \t]*complex"))
;	    (structured-types '("structure" "union" "map"))
;	    (other-types '("record" "dimension" "parameter" "common" "save"
;			   "external" "intrinsic" "data" "equivalence")))
;	(append
;	 (mapcar (lambda (x) (concat "implicit[ \t]*" x)) simple-types)
;	 simple-types
;	 (mapcar (lambda (x) (concat "end[ \t]*" x)) structured-types)
;	 structured-types
;	 other-types)))
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
	 (fortran-keywords
;	  ("continue" "format" "end" "enddo" "if" "then" "else" "endif"
;	   "elseif" "while" "inquire" "stop" "return" "include" "open"
;	   "close" "read" "write" "format" "print")
	  (concat "c\\(lose\\|ontinue\\)\\|"
		  "e\\(lse\\(\\|if\\)\\|nd\\(\\|do\\|if\\)\\)\\|format\\|"
		  "i\\(f\\|n\\(clude\\|quire\\)\\)\\|open\\|print\\|"
		  "re\\(ad\\|turn\\)\\|stop\\|then\\|w\\(hile\\|rite\\)"))
        (fortran-logicals
;       ("and" "or" "not" "lt" "le" "eq" "ge" "gt" "ne" "true" "false")
         "and\\|eq\\|false\\|g[et]\\|l[et]\\|n\\(e\\|ot\\)\\|or\\|true"))

  (setq fortran-font-lock-keywords-1
   (list
    ;;
    ;; Fontify syntactically (assuming strings cannot be quoted or span lines).
    (cons (concat "^[" comment-chars "].*") 'font-lock-comment-face)
    '(fortran-match-!-comment . font-lock-comment-face)
    (list (concat "^[^" comment-chars "\t\n]" (make-string 71 ?.) "\\(.*\\)")
	  '(1 font-lock-comment-face))
    '("'[^'\n]*'?" . font-lock-string-face)
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
     (cons (concat "\\<\\(" fortran-type-types "\\)\\>") 'font-lock-type-face)
     ;;
     ;; Fontify all builtin keywords (except logical, do and goto; see below).
     (concat "\\<\\(" fortran-keywords "\\)\\>")
     ;;
     ;; Fontify all builtin operators.
     (concat "\\.\\(" fortran-logicals "\\)\\.")
     ;;
     ;; Fontify do/goto keywords and targets, and goto tags.
     (list "\\<\\(do\\|go *to\\)\\>[ \t]*\\([0-9]+\\)?"
	   '(1 font-lock-keyword-face)
	   '(2 font-lock-reference-face nil t))
     (cons "^ *\\([0-9]+\\)" 'font-lock-reference-face))))

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
	   '(font-lock-match-c++-style-declaration-item-and-skip-to-next
	     ;; Start after any *(...) expression.
	     (and (match-beginning 15) (forward-sexp 1))
	     ;; No need to clean up.
	     nil
	     ;; Fontify as a variable name, functions are fontified elsewhere.
	     (1 font-lock-variable-name-face nil t))))
    ;;
    ;; Things extra to `fortran-font-lock-keywords-3' (must be done first).
    (list
     ;;
     ;; Fontify goto-like `err=label'/`end=label' in read/write statements.
     '(", *\\(e\\(nd\\|rr\\)\\)\\> *\\(= *\\([0-9]+\\)\\)?"
       (1 font-lock-keyword-face) (4 font-lock-reference-face nil t))
     ;;
     ;; Highlight standard continuation character and in a TAB-formatted line.
     '("^     \\([^ 0]\\)" 1 font-lock-string-face)
     '("^\t\\([1-9]\\)" 1 font-lock-string-face))
    ;;
    ;; The list `fortran-font-lock-keywords-2' less that for types (see above).
    (cdr (nthcdr (length fortran-font-lock-keywords-1)
		 fortran-font-lock-keywords-2))))
  )

(defvar fortran-font-lock-keywords fortran-font-lock-keywords-1
  "Default expressions to highlight in Fortran mode.")


(defvar fortran-mode-map () 
  "Keymap used in Fortran mode.")
(if fortran-mode-map
    ()
  (setq fortran-mode-map (make-sparse-keymap))
  (define-key fortran-mode-map ";" 'fortran-abbrev-start)
  (define-key fortran-mode-map "\C-c;" 'fortran-comment-region)
  (define-key fortran-mode-map "\e\C-a" 'beginning-of-fortran-subprogram)
  (define-key fortran-mode-map "\e\C-e" 'end-of-fortran-subprogram)
  (define-key fortran-mode-map "\e;" 'fortran-indent-comment)
  (define-key fortran-mode-map "\e\C-h" 'mark-fortran-subprogram)
  (define-key fortran-mode-map "\e\n" 'fortran-split-line)
  (define-key fortran-mode-map "\n" 'fortran-indent-new-line)
  (define-key fortran-mode-map "\e\C-q" 'fortran-indent-subprogram)
  (define-key fortran-mode-map "\C-c\C-w" 'fortran-window-create-momentarily)
  (define-key fortran-mode-map "\C-c\C-r" 'fortran-column-ruler)
  (define-key fortran-mode-map "\C-c\C-p" 'fortran-previous-statement)
  (define-key fortran-mode-map "\C-c\C-n" 'fortran-next-statement)
  (define-key fortran-mode-map "\t" 'fortran-indent-line)
  (define-key fortran-mode-map "0" 'fortran-electric-line-number)
  (define-key fortran-mode-map "1" 'fortran-electric-line-number)
  (define-key fortran-mode-map "2" 'fortran-electric-line-number)
  (define-key fortran-mode-map "3" 'fortran-electric-line-number)
  (define-key fortran-mode-map "4" 'fortran-electric-line-number)
  (define-key fortran-mode-map "5" 'fortran-electric-line-number)
  (define-key fortran-mode-map "6" 'fortran-electric-line-number)
  (define-key fortran-mode-map "7" 'fortran-electric-line-number)
  (define-key fortran-mode-map "8" 'fortran-electric-line-number)
  (define-key fortran-mode-map "9" 'fortran-electric-line-number))

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

;;;###autoload
(defun fortran-mode ()
  "Major mode for editing Fortran code.
\\[fortran-indent-line] indents the current Fortran line correctly. 
DO statements must not share a common CONTINUE.

Type ;? or ;\\[help-command] to display a list of built-in\
 abbrevs for Fortran keywords.

Key definitions:
\\{fortran-mode-map}

Variables controlling indentation style and extra features:

 comment-start
    Normally nil in Fortran mode.  If you want to use comments
    starting with `!', set this to the string \"!\".
 fortran-do-indent
    Extra indentation within do blocks.  (default 3)
 fortran-if-indent
    Extra indentation within if blocks.  (default 3)
 fortran-structure-indent
    Extra indentation within structure, union, map and interface blocks.
    (default 3)
 fortran-continuation-indent
    Extra indentation applied to continuation statements.  (default 5)
 fortran-comment-line-extra-indent
    Amount of extra indentation for text within full-line comments. (default 0)
 fortran-comment-indent-style
    nil    means don't change indentation of text in full-line comments,
    fixed  means indent that text at `fortran-comment-line-extra-indent' beyond
           the value of `fortran-minimum-statement-indent-fixed' (for fixed
           format continuation style) or `fortran-minimum-statement-indent-tab'
           (for TAB format continuation style).
    relative  means indent at `fortran-comment-line-extra-indent' beyond the
 	      indentation for a line of code.
    (default 'fixed)
 fortran-comment-indent-char
    Single-character string to be inserted instead of space for
    full-line comment indentation.  (default \" \")
 fortran-minimum-statement-indent-fixed
    Minimum indentation for Fortran statements in fixed format mode. (def.6)
 fortran-minimum-statement-indent-tab
    Minimum indentation for Fortran statements in TAB format mode. (default 9)
 fortran-line-number-indent
    Maximum indentation for line numbers.  A line number will get
    less than this much indentation if necessary to avoid reaching
    column 5.  (default 1)
 fortran-check-all-num-for-matching-do
    Non-nil causes all numbered lines to be treated as possible \"continue\"
    statements.  (default nil)
 fortran-blink-matching-if 
    Non-nil causes \\[fortran-indent-line] on an ENDIF statement to blink on
    matching IF.  Also, from an ENDDO statement, blink on matching DO [WHILE]
    statement.  (default nil)
 fortran-continuation-string
    Single-character string to be inserted in column 5 of a continuation
    line.  (default \"$\")
 fortran-comment-region
    String inserted by \\[fortran-comment-region] at start of each line in 
    region.  (default \"c$$$\")
 fortran-electric-line-number
    Non-nil causes line number digits to be moved to the correct column 
    as typed.  (default t)
 fortran-break-before-delimiters
    Non-nil causes `fortran-fill' breaks lines before delimiters.
    (default t)
 fortran-startup-message
    Set to nil to inhibit message first time Fortran mode is used.

Turning on Fortran mode calls the value of the variable `fortran-mode-hook'
with no args, if that value is non-nil."
  (interactive)
  (kill-all-local-variables)
  (if fortran-startup-message
      (message "Emacs Fortran mode %s.  Bugs to %s"
	       fortran-mode-version bug-fortran-mode))
  (setq fortran-startup-message nil)
  (setq local-abbrev-table fortran-mode-abbrev-table)
  (set-syntax-table fortran-mode-syntax-table)
  ;; Font Lock mode support.
  (make-local-variable 'font-lock-defaults)
  (setq font-lock-defaults '((fortran-font-lock-keywords
			      fortran-font-lock-keywords-1
			      fortran-font-lock-keywords-2
			      fortran-font-lock-keywords-3)
			     t t ((?/ . "$/"))))
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
  (setq fill-column 72) ; Already local?
  (use-local-map fortran-mode-map)
  (setq mode-name "Fortran")
  (setq major-mode 'fortran-mode)
;;;(make-local-variable 'fortran-tab-mode)
  (make-local-variable 'fortran-comment-line-extra-indent)
  (make-local-variable 'fortran-minimum-statement-indent-fixed)
  (make-local-variable 'fortran-minimum-statement-indent-tab)
  (make-local-variable 'fortran-column-ruler-fixed)
  (make-local-variable 'fortran-column-ruler-tab)
  (make-local-variable 'fortran-tab-mode-string)
  (setq fortran-tab-mode-string " TAB-format")
  (setq indent-tabs-mode (fortran-analyze-file-format))
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
if the value of  comment-start  is not nil.
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
		      (- (calculate-fortran-indent) (current-column))))))

(defun fortran-comment-region (beg-region end-region arg)
  "Comments every line in the region.
Puts fortran-comment-region at the beginning of every line in the region. 
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
  "Inserts a column ruler momentarily above current line, till next keystroke.
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
  "Makes the window 72 columns wide.
See also `fortran-window-create-momentarily'."
  (interactive)
  (condition-case error
      (progn
	(let ((window-min-width 2))
	  (if (< (window-width) (frame-width))
	      (enlarge-window-horizontally (- (frame-width)
					      (window-width) 1)))
	  (split-window-horizontally 73)
	  (other-window 1)
	  (switch-to-buffer " fortran-window-extra" t)
	  (select-window (previous-window))))
    (error (message "No room for Fortran window.")
	   'error)))

(defun fortran-window-create-momentarily (&optional arg)
  "Momentarily makes the window 72 columns wide.
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
	(progn 
	  (insert "\n\t")
	  (insert-char (fortran-numerical-continuation-char) 1))
      (insert "\n " fortran-continuation-string)));Space after \n important
  (fortran-indent-line))		;when the cont string is C, c or *.

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

(defun fortran-electric-line-number (arg)
  "Self insert, but if part of a Fortran line number indent it automatically.
Auto-indent does not happen if a numeric arg is used."
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
				  t)) ;not a line number
	    (looking-at "[0-9]")	;within a line number
	    )
	(self-insert-command (prefix-numeric-value arg))
      (skip-chars-backward " \t")
      (insert last-command-char)
      (fortran-indent-line))))

(defun beginning-of-fortran-subprogram ()
  "Moves point to the beginning of the current Fortran subprogram."
  (interactive)
  (let ((case-fold-search t))
    (beginning-of-line -1)
    (re-search-backward "^[ \t0-9]*end\\b[ \t]*[^ \t=(a-z]" nil 'move)
    (if (looking-at "^[ \t0-9]*end\\b[ \t]*[^ \t=(a-z]")
	(forward-line 1))))

(defun end-of-fortran-subprogram ()
  "Moves point to the end of the current Fortran subprogram."
  (interactive)
  (let ((case-fold-search t))
    (beginning-of-line 2)
    (re-search-forward "^[ \t0-9]*end\\b[ \t]*[^ \t=(a-z]" nil 'move)
    (goto-char (match-beginning 0))
    (forward-line 1)))

(defun mark-fortran-subprogram ()
  "Put mark at end of Fortran subprogram, point at beginning. 
The marks are pushed."
  (interactive)
  (end-of-fortran-subprogram)
  (push-mark (point))
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

(defun fortran-blink-matching-if ()
  ;; From a Fortran ENDIF statement, blink the matching IF statement.
  (let ((top-of-window (window-start)) matching-if
	(endif-point (point)) message)
    (if (save-excursion (beginning-of-line)
			(skip-chars-forward " \t0-9")
			(looking-at "end[ \t]*if\\b"))
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
  ;; From a Fortran ENDDO statement, blink on the matching DO or DO WHILE
  ;; statement.  This is basically copied from fortran-blink-matching-if.
  (let ((top-of-window (window-start)) matching-do
	(enddo-point (point)) message)
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
  ;; Search forward for first unmatched ENDDO.  Return point or nil.
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
                    (not (looking-at "^[ \t0-9]*end\\b[ \t]*[^ \t=(a-z]")))

          (skip-chars-forward " \t0-9")
          (cond ((looking-at "end[ \t]*do\\b")
                 (setq count (- count 1)))
                ((looking-at "do[ \t]+[^0-9]")
                 (setq count (+ count 1)))))
        (and (= count 0)
             ;; All pairs accounted for.
             (point))))))

(defun fortran-beginning-do ()
  ;; Search backwards for first unmatched DO [WHILE].  Return point or nil.
  (if (save-excursion (beginning-of-line)
                      (skip-chars-forward " \t0-9")
                      (looking-at "do[ \t]+"))
      ;; Sitting on one.
      (match-beginning 0)
    ;; Search for one.
    (save-excursion
      (let ((count 1))
        (while (and (not (= count 0))
                    (not (eq (fortran-previous-statement) 'first-statement))
                    ;; Keep local to subprogram
                    (not (looking-at "^[ \t0-9]*end\\b[ \t]*[^ \t=(a-z]")))

          (skip-chars-forward " \t0-9")
          (cond ((looking-at "do[ \t]+[^0-9]")
                 (setq count (- count 1)))
                ((looking-at "end[ \t]*do\\b")
                 (setq count (+ count 1)))))

        (and (= count 0)
             ;; All pairs accounted for.
             (point))))))

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

(defun fortran-end-if ()
  ;; Search forwards for first unmatched ENDIF.  Return point or nil.
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
                    (not (looking-at
                          "^[ \t0-9]*end\\b[ \t]*[^ \t=(a-z]")))

          (skip-chars-forward " \t0-9")
          (cond ((looking-at "end[ \t]*if\\b")
                 (setq count (- count 1)))

                ((looking-at "if[ \t]*(")
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
             (point))))))

(defun fortran-beginning-if ()
  ;; Search backwards for first unmatched IF-THEN.  Return point or nil.
  (if (save-excursion
        ;; May be sitting on multi-line if-then statement, first move to
        ;; beginning of current statement.  Note: `fortran-previous-statement'
        ;; moves to previous statement *unless* current statement is first
        ;; one.  Only move forward if not first-statement.
        (if (not (eq (fortran-previous-statement) 'first-statement))
            (fortran-next-statement))
        (skip-chars-forward " \t0-9")
        (and
         (looking-at "if[ \t]*(")
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
                    (not (looking-at
                          "^[ \t0-9]*end\\b[ \t]*[^ \t=(a-z]")))

          (skip-chars-forward " \t0-9")
          (cond ((looking-at "if[ \t]*(")
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
             (point))))))

(defun fortran-indent-line ()
  "Indents current Fortran line based on its contents and on previous lines."
  (interactive)
  (let ((cfi (calculate-fortran-indent)))
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
An abbrev before point is expanded if `abbrev-mode' is non-nil."
  (interactive)
  (if abbrev-mode (expand-abbrev))
  (save-excursion
    (beginning-of-line)
    (skip-chars-forward " \t")
    (if (or (looking-at "[0-9]")	;Reindent only where it is most
	    (looking-at "end")		;likely to be necessary
	    (looking-at "else")
	    (looking-at (regexp-quote fortran-continuation-string)))
	(fortran-indent-line)))
  (newline)
  (fortran-indent-line))

(defun fortran-indent-subprogram ()
  "Properly indents the Fortran subprogram which contains point."
  (interactive)
  (save-excursion
    (mark-fortran-subprogram)
    (message "Indenting subprogram...")
    (indent-region (point) (mark) nil))
  (message "Indenting subprogram...done."))

(defun calculate-fortran-indent ()
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
	  (cond ((looking-at "if[ \t]*(")
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
		((looking-at "\\(else\\|elseif\\)\\b")
		 (setq icol (+ icol fortran-if-indent)))
		((looking-at "select[ \t]*case[ \t](.*)\\b")
		 (setq icol (+ icol fortran-if-indent)))
		((looking-at "case[ \t]*(.*)[ \t]*\n")
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
		((looking-at "end\\b[ \t]*[^ \t=(a-z]")
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
	     (cond ((looking-at "end[ \t]*if\\b")
		    (setq icol (- icol fortran-if-indent)))
		   ((looking-at "\\(else\\|elseif\\)\\b")
		    (setq icol (- icol fortran-if-indent)))
                   ((looking-at "case[ \t]*(.*)[ \t]*\n")
		    (setq icol (- icol fortran-if-indent)))
                   ((looking-at "case[ \t]*default\\b")
		    (setq icol (- icol fortran-if-indent)))
		   ((looking-at "\\(otherwise\\|else[ \t]*where\\)\\b")
		    (setq icol (- icol fortran-if-indent)))
		   ((looking-at "end[ \t]*where\\b")
		    (setq icol (- icol fortran-if-indent)))
		   ((and (looking-at "continue\\b")
			 (fortran-check-for-matching-do))
		    (setq icol (- icol fortran-do-indent)))
		   ((looking-at "end[ \t]*do\\b")
		    (setq icol (- icol fortran-do-indent)))
		   ((looking-at
		     "end[ \t]*\
\\(structure\\|union\\|map\\|interface\\)\\b[ \t]*[^ \t=(a-z]")
		    (setq icol (- icol fortran-structure-indent)))
		   ((looking-at
		     "end[ \t]*select\\b[ \t]*[^ \t=(a-z]")
		    (setq icol (- icol fortran-if-indent)))
		   ((and (looking-at "end\\b[ \t]*[^ \t=(a-z]")
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
  "Indents current line with spaces to column COL.
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
  "When called from a numbered statement, returns t if matching DO is found.
Otherwise return a nil."
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
		 (looking-at (concat "^[ \t0-9]*do[ \t]*0*" charnum))))))))

(defun fortran-find-comment-start-skip ()
  "Move to past `comment-start-skip' found on current line.
Return t if `comment-start-skip' found, nil if not."
;;; In order to move point only if comment-start-skip is found,
;;; this one uses a lot of save-excursions.  Note that re-search-forward
;;; moves point even if comment-start-skip is inside a string-constant.
;;; Some code expects certain values for match-beginning and end
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

;;;From: simon@gnu (Simon Marshall)
;;; Find the next ! not in a string.
(defun fortran-match-!-comment (limit)
  (let (found)
    (while (and (setq found (search-forward "!" limit t))
                (fortran-is-in-string-p (point))))
    (if (not found)
	nil
      ;; Cheaper than `looking-at' "!.*".
      (store-match-data
       (list (1- (point)) (progn (end-of-line) (min (point) limit))))
      t)))

;; The above function is about 10% faster than the below...
;;(defun fortran-match-!-comment (limit)
;;  (let (found)
;;    (while (and (setq found (re-search-forward "!.*" limit t))
;;                (fortran-is-in-string-p (match-beginning 0))))
;;    found))

;;;From: ralf@up3aud1.gwdg.de (Ralf Fassel)
;;; Test if TAB format continuation lines work.
(defun fortran-is-in-string-p (where)
  "Return non-nil if POS (a buffer position) is inside a Fortran string,
nil else."
  (save-excursion
    (goto-char where)
    (cond
     ((bolp) nil)			; bol is never inside a string
     ((save-excursion			; comment lines too
	(beginning-of-line)(looking-at comment-line-start-skip)) nil)
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
		;; get around a bug in forward-line in versions <= 18.57
		(if (or (> (forward-line 1) 0) (eobp))
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
		   'fortran-do-auto-fill
		 nil))
    (force-mode-line-update)))

(defun fortran-do-auto-fill ()
  (if (> (current-column) fill-column)
      (fortran-indent-line)))

(defun fortran-fill ()
  (interactive)
  (let* ((opoint (point))
	 (bol (save-excursion (beginning-of-line) (point)))
	 (eol (save-excursion (end-of-line) (point)))
	 (bos (min eol (+ bol (fortran-current-line-indentation))))
	 (quote
	  (save-excursion
	    (goto-char bol)
	    (if (looking-at comment-line-start-skip)
		nil			; OK to break quotes on comment lines.
	      (move-to-column fill-column)
	      (cond ((fortran-is-in-string-p (point))
		     (save-excursion (re-search-backward "[^']'[^']" bol t)
				     (if fortran-break-before-delimiters
					 (point)
				       (1+ (point)))))
		    (t nil)))))
	 ;;
	 ;; decide where to split the line. If a position for a quoted
	 ;; string was found above then use that, else break the line
	 ;; before the last delimiter.
	 ;; Delimiters are whitespace, commas, and operators.
	 ;; Will break before a pair of *'s.
	 ;;
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
;;;what is this doing???
		      (if (not (re-search-forward "[\t\n,'+-/*)=]" eol t))
			  (goto-char bol))))
		(if (bolp)
		    (re-search-forward "[ \t]" opoint t)
		  (forward-char -1)
		  (if (looking-at "'")
		      (forward-char 1)
		    (skip-chars-backward " \t\*")))
		(if fortran-break-before-delimiters
		    (point)
		  (1+ (point))))))
	 )
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
		   (+ (calculate-fortran-indent) fortran-continuation-indent))
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
;;; Forward line 1 really needs to go to next non white line
    (if (save-excursion (forward-line 1)
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
  "Returns nil if fixed format is used, t if TAB formatting is used.
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

; --- The following are extensions by Mark van Schilfgaarde ---
(defun fortran-first-executable-statement ()
"  Moves point to first executable statement in current subprogram."
   (interactive)
   (end-of-fortran-subprogram)
   (beginning-of-fortran-subprogram) (fortran-next-statement)
   (while (looking-at
            " +\\(character\\|common\\|complex\\|data\\|dimension\\|double\\|equivalence\\|external\\|function\\|implicit\\|integer\\|intrinsic\\|logical\\|parameter\\|program\\|real\\|save\\|subroutine\\)\\b")
      (fortran-next-statement))
)


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
     (if (looking-at "^\\. ") (progn (delete-char 2) (insert "\n")))
     (end-of-line) (while (< (current-column) sizcol) (insert " "))
     (setq thispt (point)) (insert subrstrng) (setq endpt (point))
     (goto-char thispt) (forward-sexp)

     (while (looking-at " *[(a-z0-9+*/]") (forward-sexp))
; ... Skip past floating point numbers
     (while (looking-at "[0-9.EeDd-+]") (forward-char 1))
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

(defun fortran-next-block ()
"  Executes fortran-next-statement until fortran-end-of-block is t.
   However, if on a 'do' statement that has a matching continuation,
   go to the matching continuation."

   (interactive)

; Cursor at beginning of current statement
   (fortran-next-statement) (fortran-previous-statement)
   (setq first-statement nil)

   (while
   (cond

; If at beginning of do loop, go to the end
     ((looking-at "......[ \t]*do[ \t]+[0-9]")
      (if (not first-statement)
      (progn
	(re-search-forward "......[ \t]*do[ \t]+0*")
	(setq charnum (buffer-substring (point) (progn (skip-chars-forward "0-9") (point))))
	(if (not (re-search-forward
		   (concat "\\(^[ \t0-9]*end\\b[ \t]*$\\)\\|\\(^[ \t]*"
			   charnum "\\b\\)")
		   nil t)) (goto-char (point-max)))
	(beginning-of-line))) nil)

; If at beginning of an "if-then" block ...
     ((and (looking-at "......[ \t]*if[ \t]*(")
	   (save-excursion (forward-sexp 2)
		  (looking-at "[ \t]*then[ \t]*$")))
      (if (not first-statement) (progn
      (fortran-next-statement) (fortran-end-of-if-block t))) nil)

     ((looking-at "......[ \t]*else")
      (if (not first-statement) (progn
      (fortran-next-statement) (fortran-end-of-if-block t))) nil)

     ((fortran-check-for-matching-do)
      (if (not first-statement)
	  (progn (setq first-statement t) (fortran-next-statement) t)
	  nil))

     ((progn (setq first-statement t) (fortran-next-statement) t))
   ))
)

(defun fortran-previous-block ()
"  Executes fortran-previous-statement until fortran-end-of-block is t.
   However, if on a 'do' statement that has a matching continuation,
   go to the matching continuation."

   (interactive)

; Cursor at beginning of current statement
   (fortran-next-statement) (fortran-previous-statement)
   (setq first-statement nil)

   (while
   (cond

     ((fortran-check-for-matching-do)
      (if (not first-statement)
       (progn
	 (skip-chars-forward " \t")
	 (skip-chars-forward "0") ;skip past leading zeros
	 (setq charnum
	       (buffer-substring (point)
				 (progn (skip-chars-forward "0-9") (point))))
	 (re-search-backward (concat "^[ \t0-9]*do[ \t]*0*" charnum))
	 (beginning-of-line)
	 )) nil)

; If at beginning bottom of an "if-then" block ...
     ((looking-at "......[ \t]*end[ \t]*if")
      (if (not first-statement) (progn
      (fortran-previous-statement) (fortran-start-of-if-block t))) nil)

     ((looking-at "......[ \t]*else")
      (if (not first-statement) (progn
      (fortran-previous-statement) (fortran-start-of-if-block t))) nil)

; If at beginning of do loop ...
     ((looking-at "......[ \t]*do[ \t]+[0-9]")
      (if (not first-statement)
	  (progn (setq first-statement t) (fortran-previous-statement) t)
	  nil))

     ((progn (setq first-statement t) (fortran-previous-statement) t))
   ))
)

(defun fortran-end-of-if-block (ifelse)
  (interactive)

; Cursor at beginning of current statement
   (fortran-next-statement) (fortran-previous-statement)

   (while (and (setq not-last-statement (save-excursion (= (forward-line 1) 0)))
	       (not
		 (or (looking-at "......[ \t]*end[ \t]*if")
		     (and ifelse (looking-at "......[ \t]*else"))
		     (looking-at "......[ \t]*end[ \t]*$"))))
     (if (and (looking-at "......[ \t]*if[ \t]*(")
	      (save-excursion (forward-sexp 2)
			      (looking-at "[ \t]*then[ \t]*$")))
	 (progn (fortran-next-statement)
		(fortran-end-of-if-block nil)))
	 (fortran-next-statement))

   t
)

(defun fortran-start-of-if-block (ifelse)
  (interactive)

; Cursor at beginning of current statement
   (fortran-next-statement) (fortran-previous-statement)

   (while (and (save-excursion (= (forward-line -1) 0))
	       (not
		 (or (and (looking-at "......[ \t]*if[ \t]*(")
			  (save-excursion
			    (forward-sexp 2)
			    (looking-at "[ \t]*then[ \t]*$")))
		     (and ifelse (looking-at "......[ \t]*else"))
		     (looking-at "......[ \t]*subroutine"))))
     (if (looking-at "......[ \t]*end[ \t]*if")
	 (progn (fortran-previous-statement)
		(fortran-start-of-if-block nil)))
	 (fortran-previous-statement))
   t
)
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

  " Called by fortran-call-tree"

  (interactive "P")

  (let (tag-already-found tagnam strn)
; --- For each call do ---
  (while (re-search-forward "^[ 0-9]....+\\bcall +" eos t)
    (progn
    (setq pointnow (point))
    (setq tagnam
	  (buffer-substring
	    (point) (save-excursion
		      (skip-chars-forward "A-Za-z0-9") (point))))
    (setq tagfil (if (is-tag tagnam)
		     (progn 
		       (visit-tags-table-buffer)
		       (file-of-tag))
		      "*"))
    (set-buffer "*Tree*") (goto-char (point-max))
    (setq tag-already-found
	  (save-excursion
	    (re-search-backward
	      (concat "\\b" tagnam "\\b") (point-min) t)))
;   (snit)
    (setq strn (concat tagnam "  "
		       (if tag-already-found " ... " tagfil)))
    (insert spaces strn "\n")
    (set-buffer bnam) (setq pointnow (point))
;    (message "%s%s" spaces strn)
    (if (not outmode) (message "%s%s" spaces strn))
    )

    (if (and (is-tag tagnam) (> level 0) (not tag-already-found))
	(progn
	  (find-tagx tagnam)
          (if outmode (message "%s" (buffer-file-name)))
	  (setq buffer (car (buffer-list)))
	  (setq tnam (buffer-name buffer))
	  (fortran-call-tree2 tnam (- level 1) outmode
               	 (concat "  " spaces) (point)
	     (save-excursion (end-of-fortran-subprogram) (point)))
	  (pop-to-buffer bnam) (goto-char pointnow)))))
)

(defun find-tagx (tagname &optional next other-window)
  "Indentical to function FIND-TAG, except that find-tagx searches
   the regular expression \\bTAGNAME\\b in the tags table,
   which avoids accidental tags with longer names."

  (interactive (if current-prefix-arg
		   '(nil t)
		 (find-tag-tag "Find tag: ")))
  (let (buffer file linebeg startpos)
    (save-excursion
     (visit-tags-table-buffer)
     (if (not next)
	 (goto-char (point-min))
       (setq tagname last-tag))
     (setq last-tag tagname)
     (while (progn
	      (if (not (re-search-forward (concat "\\b" tagname "\\b") nil t))
		  (error "No %sentries containing %s"
			 (if next "more " "") tagname))
	      (not (looking-at "[^\n\177]*\177"))))
     (search-forward "\177")
     (setq file (expand-file-name (file-of-tag)
				  (file-name-directory tags-file-name)))
     (setq linebeg
	   (buffer-substring (1- (point))
			     (save-excursion (beginning-of-line) (point))))
     (search-forward ",")
     (setq startpos (read (current-buffer))))
    (if other-window
	(find-file-other-window file)
      (find-file file))
    (widen)
    (push-mark)
    (let ((offset 1000)
	  found
	  (pat (concat "^" (regexp-quote linebeg))))
      (or startpos (setq startpos (point-min)))
      (while (and (not found)
		  (progn
		   (goto-char (- startpos offset))
		   (not (bobp))))
	(setq found
	      (re-search-forward pat (+ startpos offset) t))
	(setq offset (* 3 offset)))
      (or found
	  (re-search-forward pat nil t)
	  (error "%s not found in %s" pat file)))
    (beginning-of-line))
  (setq tags-loop-form '(find-tag nil t))
  ;; Return t in case used as the tags-loop-form.
  t)

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
"  Forward past any f77 'white space', which includes blanks, tabs,
   comment lines and the first 6 columns in a continuation line.
   Optional ARG non-nil includes characters -+*/=, as whitespace.
   Returns t if position point changes, otherwise nil."

   (interactive "P")

   (let ((pointnow (point)) (wschar "-+*/=, \t\n"))
     (if (not arg) (setq wschar " \t\n"))
     (while
	 (or (= (current-column) 5)
	     (looking-at (concat "[" wschar "]"))
	     (fortran-is-comment-line))
					; Skip over comments
       (if (fortran-is-comment-line) (fortran-next-statement))
					; and over whitespace
       (skip-chars-forward wschar)
					; and over the continuation char
       (if (= (current-column) 5) (forward-char 1)))
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
       (if (fortran-is-comment-line)
	   (progn (fortran-previous-statement) (forward-line 1)))
					; and over whitespace
       (skip-chars-backward wschar)
					; and over the continuation char
       (if (= (current-column) 6) (backward-char 1)))
     (not (= (point) pointnow)))
)

; --- macro fortran-forward-variable ---
(defun fortran-forward-variable (&optional arg)
"  Forward past one variable,array or function.
   Invokes fortran-forward-whitespace; then moves past a word
   and subsequently past a sexp if it follows point.
   Optional ARG repeats ARG times."

   (interactive "p") (if (not arg) (setq arg 1))

   (let ((vchar "A-Za-z0-9._"))
     (while (> arg 0)
       (progn
	 (fortran-forward-whitespace 1)
	 (if (looking-at ")")
	     (progn (up-list 1) (fortran-forward-whitespace 1)))
	 (progn
	   (skip-chars-forward vchar)
;(while (looking-at (concat "[ \t]*[" vchar "]"))
;(progn
;(fortran-forward-whitespace 0) (skip-chars-forward vchar)))
	   (if (looking-at "[ \t]*(") (forward-sexp))
	   (if (looking-at "'") (forward-sexp)))
	 (setq arg (- arg 1)))))
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

; --- macro fortran-is-comment-line ---
(defun fortran-is-comment-line (&optional arg)
"  Returns t if point in a fortran comment, else nil"

   (interactive "p")
   (save-excursion (beginning-of-line) (looking-at "^[Cc*]"))
)

(defun compilation-sgi-declare-fortran-variable-type ()
"  Get next undeclared variable in compilation buffer and insert
   a declaration in prior window"
   (interactive)

   (re-search-forward "Variable \"\\(\\w+\\)\" must be explicitly typed")
   (setq vnam (buffer-substring (match-beginning 1) (match-end 1)))
   (setq bufnow  (buffer-name)
         pointnow (point)
         buffer (find-file-noselect (car gud-last-last-frame))
 	 window (display-buffer bufnow))

   (other-window -1) (fortran-first-executable-statement)
   (if (and (>= (string-to-char vnam) ?I) (<= (string-to-char vnam) ?O))
       (re-search-backward "^ +integer")
       (re-search-backward "^ +double precision"))
   (fortran-next-statement) (backward-char 1)
   (insert-string (concat "," vnam))
   (backward-word 1) (downcase-word 1)
   (other-window 1)
   (message "inserting variable %s" vnam)

)

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
         (progn (backward-word 1) (delete-char 5) (insert-string "min")
              (setq posn (- posn 2))
   ))))

; Handle nint,min,max
    (goto-char startposn)
    (while (re-search-forward
            "\\b\\(nint\\|min\\|max\\)\\b"  posn t) (progn
      (setq query (get-query "Convert?" query))
      (if (or (char-equal query ?a) (char-equal query ?y))
         (progn (backward-word 1) (insert-string "x")
              (setq posn (+ posn 1))
   ))))

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

(provide 'fortran)

;;; fortran.el ends here

