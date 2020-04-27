(defun ctrlv6to7 ()
" Approximately convert ctrl file v6 to v7"
    (interactive)

  (let ((pointnow (point))
	strn query f)

    ;; Add to VERS
    (goto-char (point-min)) (setq query ?y) (setq query (get-query "Add to VERS?" query))
    (if (or (char-equal query ?a) (char-equal query ?y))
	(if (re-search-forward "^VERS " (point-max) f)
	    (progn (end-of-line) (insert-string " LM:7 ASA:7")
		   (setq query ?y) (setq query (get-query "edit?" query))
		   (if (or (char-equal query ?a) (char-equal query ?y))
		       (recursive-edit)))))

    ;; change NCLASS to NSPEC
    (if (save-excursion (goto-char (point-min)) (re-search-forward "NCLASS=" (point-max) t))
	(progn (setq query ?y) (setq query (get-query "Replace NCLASS with NSPEC?" query))
	       (if (or (char-equal query ?a) (char-equal query ?y))
		   (query-replace "NCLASS="  "NSPEC="))))

    ;; replace MIX with ITER
    (goto-char (point-min)) (setq query ?y) (setq query (get-query "Replace MIX with ITER?" query))
    (if (or (char-equal query ?a) (char-equal query ?y))
    (while (re-search-forward "^MIX" (point-max) t)
      (progn (beginning-of-line) (kill-line 1) (yank) (yank) (previous-line 2)
	     (delete-char 4) (insert-string "ITER")
	     (if (re-search-forward "MODE=" (save-excursion (re-search-forward "^[A-Za-z]" (point-max) t)) t)
		 (replace-match "MIX=" nil nil))
	     (if
		 (save-excursion
		   (and (re-search-forward "^START " (point-max) f)
			(re-search-forward " NIT=[^ ]+" (save-excursion (re-search-forward "^[A-Za-z]" (point-max) t)) t)))
		 (progn (end-of-line) (insert-string (match-string 0))))

	     (setq query ?y) (setq query (get-query "edit?" query))
	     (if (or (char-equal query ?a) (char-equal query ?y)) (recursive-edit))
	     )
      (re-search-forward "^MIX" (point-max) t)))

    ;; copy START_CNVG to ITER_CONVC
    (if
	(save-excursion
	  (goto-char (point-min))
	  (and (re-search-forward "^ITER " (point-max) t) 
	       (setq pointnow (point))
	       (re-search-forward "^START " (point-max) t)
	       (re-search-forward " CNVG=\\([^ ]+\\)" (save-excursion (re-search-forward "^[A-Za-z]" (point-max) t)) t)))
	(progn (setq query ?y) (setq query (get-query "Copy START_CNVG to ITER_CONVC?" query))
	       (if (or (char-equal query ?a) (char-equal query ?y))
		   (progn (goto-char pointnow)
			  (end-of-line)
			  (insert-string (concat " CONVC=" (match-string 1)))))))

    ;; New category HAM
    (if
	(save-excursion
	  (goto-char (point-min))
	  (not (re-search-forward "^HAM " (point-max) t)))
	(progn (setq query ?y) (setq query (get-query "Create new HAM category?" query))
	       (if (or (char-equal query ?a) (char-equal query ?y))
		   (progn
		     (goto-char (point-min))  ; put HAM just before OPTIONS; else at start of file
		     (if (re-search-forward "^OPTIONS " (point-max) t)
			 (beginning-of-line)
		       (goto-char (point-min)))
		     (insert-string "HAM     \n")))))

		     
    ;; copy OPTIONS_NSPIN to HAM_NSPIN
    (if
	(save-excursion
	  (goto-char (point-min))
	  (and (re-search-forward "^HAM\\b" (point-max) t)
	       (setq pointnow (point))
	       (not (re-search-forward " NSPIN=\\([^ ]+\\)" (save-excursion (re-search-forward "^[A-Za-z]" (point-max) t)) t))
	       (goto-char (point-min))
	       (re-search-forward "^OPTIONS\\b" (point-max) t)
	       (re-search-forward " NSPIN=\\([^ ]+\\)" (save-excursion (re-search-forward "^[A-Za-z]" (point-max) t)) t)))
	(progn (setq query ?y) (setq query (get-query "Copy OPTIONS_NSPIN to HAM_NSPIN?" query))
	       (if (or (char-equal query ?a) (char-equal query ?y))
		   (progn (goto-char pointnow)
			  (end-of-line)
			  (insert-string (concat " NSPIN=" (match-string 1)))))))

    ;; copy OPTIONS_REL to HAM_REL
    (if
	(save-excursion
	  (goto-char (point-min))
	  (and (re-search-forward "^HAM\\b" (point-max) f)
	       (setq pointnow (point))
	       (not (re-search-forward " REL=\\([^ ]+\\)" (save-excursion (re-search-forward "^[A-Za-z]" (point-max) f)) t))
	       (goto-char (point-min))
	       (re-search-forward "^OPTIONS\\b" (point-max) f)
	       (re-search-forward " REL=\\([^ ]+\\)" (save-excursion (re-search-forward "^[A-Za-z]" (point-max) t)) t)))
	(progn (setq query ?y) (setq query (get-query "Copy OPTIONS_REL to HAM_REL?" query))
	       (if (or (char-equal query ?a) (char-equal query ?y))
		   (progn (goto-char pointnow)
			  (end-of-line)
			  (insert-string (concat " REL=" (match-string 1)))))))

    ;; copy OPTIONS_NONCOL to HAM_NONCOL
    (if
	(save-excursion
	  (goto-char (point-min))
	  (and (re-search-forward "^HAM\\b" (point-max) t) 
	       (setq pointnow (point))
	       (not (re-search-forward " NONCOL=\\([^ ]+\\)" (save-excursion (re-search-forward "^[A-Za-z]" (point-max) t)) t))
	       (goto-char (point-min))
	       (re-search-forward "^OPTIONS\\b" (point-max) t)
	       (re-search-forward " NONCOL=\\([^ ]+\\)" (save-excursion (re-search-forward "^[A-Za-z]" (point-max) t)) t)))
	(progn (setq query ?y) (setq query (get-query "Copy OPTIONS_NONCOL to HAM_NONCOL?" query))
	       (if (or (char-equal query ?a) (char-equal query ?y))
		   (progn (goto-char pointnow)
			  (end-of-line)
			  (insert-string (concat " NONCOL=" (match-string 1)))))))

    ;; Add new HAM QASA
    (if
	(save-excursion
	  (goto-char (point-min))
	  (and (re-search-forward "^HAM\\b" (point-max) t)
	       (setq pointnow (point))
	       (not (re-search-forward " QASA=" (save-excursion (re-search-forward "^[A-Za-z]" (point-max) t)) t))))
	(progn (setq query ?y) (setq query (get-query "Add HAM QASA=0?" query))
	       (if (or (char-equal query ?a) (char-equal query ?y))
		   (progn (goto-char pointnow)
			  (end-of-line)
			  (insert-string " QASA=0")))))

    ;; split BZ N.W
;    (snit)
    (if
	(save-excursion
	  (goto-char (point-min))
	  (and (re-search-forward "^BZ\\b" (point-max) t)
	       (setq pointnow (point))
	       (not (re-search-forward " W=" (save-excursion (re-search-forward "^[A-Za-z]" (point-max) t)) t))
	       (goto-char pointnow)
	       (re-search-forward " N[.]W= *\\([0-9]\\)[.]\\([0-9]+\\)" (save-excursion (re-search-forward "^[A-Za-z]" (point-max) t)) t)))
	(progn (setq query ?y) (setq query (get-query "Split N.W?" query))
	       (if (or (char-equal query ?a) (char-equal query ?y))
		   (progn (replace-match " N=\\1 W=0.\\2" nil nil)))))

))
