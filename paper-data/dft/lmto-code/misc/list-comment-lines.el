; --- macro list-comment-lines (25 Nov 94) ---
; Put in your .emacs, for example
; (setq fortran-mode-hook
;       '(lambda ()
; 	(make-local-variable 'fortran-mode) (setq fortran-mode t)
; 	(setq comment-line-regex2 "^[Cc] *\\(--+\\)\\|\\(\\.\\.+\\.+\\)[ ()a-zA-Z0-9]\\|^ +subroutine")
; 	(setq comment-line-regex "^[Cc] *--+[ ()a-zA-Z0-9]\\|^ +subroutine")))

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
