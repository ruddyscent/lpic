(defvar lpic-license
"
   This file is part of LPIC++, a particle-in-cell code for 
   simulating the interaction of laser light with plasma.

   Copyright (C) 1994-1997 Roland Lichters

   LPIC++ is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License
   as published by the Free Software Foundation; either version 2
   of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
")

; skeleton for a new header file---you might want to bind this to some key
(defun lpic-new-header ()
  (interactive)
  (lpic-add-license)
  (lpic-add-file-description)
  (insert "\n")
  (lpic-add-include-guard)
  (lpic-add-sample-header)
  (insert "\n")
  (lpic-add-namespace))

; skeleton for a new source file---you might want to bind this to some key
(defun lpic-new-source ()
  (interactive)
  (lpic-add-license)
  (insert "\n")
  (lpic-add-sample-header)
  (insert "\n")
  (lpic-add-namespace))

(require 'cc-mode)
(defvar c++-font-lock-extra-types)

; A few types for syntax-highlighting
; from QuantLib:
(setq c++-font-lock-extra-types
      (append c++-font-lock-extra-types
              '("QuantLib"
                "Integer" "BigInteger" "Natural" "BigNatural" "Real" "Decimal"
                "Time" "Rate" "Spread" "DiscountFactor" "Size" "Volatility"
                "Date" "Day" "Month" "Year" "Weekday"
                "TimeUnit" "Frequency" "Compounding" "Period" "DayCounter"
                "Calendar" "BusinessDayConvention"
                "Currency" "ExchangeRate" "Money" "Rounding"
                "InterestRate"
                "History"
                "Handle")))

; helper functions

(defun lpic-add-license ()
  (let ((holder (read-from-minibuffer "Copyright holder? ")))
    (let ((copyright-notice
           (apply 'string (append "   Copyright (C) "
                                  (substring (current-time-string) -4)
                                  " "
                                  holder
                                  "\n"
                                  ()))))
      (insert "\n/*\n"
              copyright-notice
              lpic-license
              "*/\n\n"))))

(defun lpic-add-file-description ()
  (let ((filename (buffer-name))
        (description (read-from-minibuffer "Short file description? ")))
    (insert "/*! \\file " filename "\n"
            "    \\brief " description "\n"
            "*/\n")))

(defun lpic-add-include-guard ()
  (let ((guard (read-from-minibuffer "Include guard? ")))
    (insert "#ifndef " guard "\n"
            "#define " guard "\n"
            "\n\n\n"
            "#endif\n"))
  (previous-line 3))

(defun lpic-add-sample-header ()
  (insert "#include <common.h>\n"))


(defun lpic-add-namespace ()
  (insert "namespace LPIC {\n"
          "\n\n\n"
          "}\n")
  (previous-line 3)
  (c-indent-command))

