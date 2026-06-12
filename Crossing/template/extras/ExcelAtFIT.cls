% ---------------------------------------------------------------------
% Excel@FIT latex template
% ---------------------------------------------------------------------
%  HISTORY:
%    2015-02-04 - Herout: intial version for Excel@FIT 2015
%    2015-04-03 - Herout: works with [czech] and [slovak] options
% ---------------------------------------------------------------------

\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{ExcelAtFIT}[25/01/2012, v1.0]
\RequirePackage{ifthen}
\RequirePackage{calc}
\AtEndOfClass{\RequirePackage{microtype}}
\DeclareOption*{\PassOptionsToClass{\CurrentOption}{article}}
\ProcessOptions*
\LoadClass[fleqn,11pt]{article}
\RequirePackage{ifpdf}      % Needed to pick between latex and pdflatex
\emergencystretch=1em

%----------------------------------------------------------------------
%	FONTS
%----------------------------------------------------------------------

\RequirePackage[T1]{fontenc}
\RequirePackage{sourcesanspro}
\renewcommand{\familydefault}{\sfdefault}

%----------------------------------------------------------------------
%	VARIOUS USEFUL PACKAGES
%----------------------------------------------------------------------

\RequirePackage[utf8]{inputenc}
\RequirePackage{amsmath,amsfonts,amssymb}
\RequirePackage{graphicx,xcolor}
\RequirePackage[english]{babel}
\RequirePackage{booktabs}
\RequirePackage{listings}

%----------------------------------------------------------------------
%	MARGINS	
%----------------------------------------------------------------------
	  	  
\RequirePackage[left=1.7cm,%
                right=1.7cm,%
				top=2.05cm,%
				bottom=2.05cm,%
				headheight=11pt,%
				a4paper]{geometry}%
				
%----------------------------------------------------------------------
%	FIGURES AND TABLES CAPTIONS
%----------------------------------------------------------------------

\RequirePackage[labelfont={bf,sf,small},%
                labelsep=period,%
                justification=raggedright]{caption}
\setlength{\abovecaptionskip}{4pt}
\setlength{\belowcaptionskip}{-10pt}
				
%----------------------------------------------------------------------
%	PAGE HEADER
%----------------------------------------------------------------------

\RequirePackage{fancyhdr}  % Needed to define custom headers/footers
\RequirePackage{lastpage}  % Number of pages in the document
\pagestyle{fancy}          % Enables the custom headers/footers
% Headers
\lhead{}%
\chead{}%
\rhead{\small\sffamily\bfseries\@PaperTitle\  --- \thepage/\pageref{LastPage}}
% Footers
\lfoot{}%
\cfoot{}%
\rfoot{}%
\renewcommand{\headrulewidth}{0pt}% % No header rule
\renewcommand{\footrulewidth}{0pt}% % No footer rule

%----------------------------------------------------------------------
%	SECTION/SUBSECTION/PARAGRAPH SET-UP
%----------------------------------------------------------------------

\RequirePackage[explicit]{titlesec}
\titleformat{\section}
  {\color{white}\large\sffamily\bfseries}
  {}
  {0em}
  {\colorbox{color-blue!100}{\parbox{\dimexpr\linewidth-2\fboxsep\relax}{\arabic{section}. #1}}}
  []
\titleformat{name=\section,numberless}
  {\color{white}\large\sffamily\bfseries}
  {}
  {0em}
  {\colorbox{color-blue!100}{\parbox{\dimexpr\linewidth-2\fboxsep\relax}{#1}}}
  []  
\titleformat{\subsection}
  {\color{color-sections}\sffamily\bfseries}
  {\thesubsection}
  {0.5em}
  {#1}
  []
\titleformat{\subsubsection}
  {\sffamily\small\bfseries}
  {\thesubsubsection}
  {0.5em}
  {#1}
  []    
\titleformat{\paragraph}[runin]
  {\sffamily\small\bfseries}
  {}
  {0em}
  {#1} 
\titlespacing*{\section}{0pc}{3ex \@plus4pt \@minus3pt}{5pt}
\titlespacing*{\subsection}{0pc}{2.5ex \@plus3pt \@minus2pt}{1pt}
\titlespacing*{\subsubsection}{0pc}{2ex \@plus2.5pt \@minus1.5pt}{0pt}
\titlespacing*{\paragraph}{0pc}{1.5ex \@plus2pt \@minus1pt}{10pt}

%----------------------------------------------------------------------
%	TABLEOFCONTENTS SET-UP
%----------------------------------------------------------------------
\newlength{\tocsep} 
\setlength\tocsep{1.5pc} % Sets the indentation of the sections in the table of contents
\setcounter{tocdepth}{3} % Three levels in the table of contents section: sections, subsections and subsubsections

\usepackage{titletoc}
\contentsmargin{0cm}
\titlecontents{section}[\tocsep]
  {\addvspace{4pt}\small\bfseries\sffamily}
  {\contentslabel[\thecontentslabel]{\tocsep}}
  {}
  {\hfill\thecontentspage}
  []
\titlecontents{subsection}[\tocsep]
  {\addvspace{2pt}\small\sffamily}
  {\contentslabel[\thecontentslabel]{\tocsep}}
  {}
  {\ \titlerule*[.5pc]{.}\ \thecontentspage}
  []
\titlecontents*{subsubsection}[\tocsep]
  {\footnotesize\sffamily}
  {}
  {}
  {}
  [\ \textbullet\ ]  

%----------------------------------------------------------------------
%	MULTIPLE AUTHOR SET
%----------------------------------------------------------------------  
  
\newcount\@authcnt
\newcount\@tmpcnt\@tmpcnt\z@

\def\@affiliation{%
  \ifnum\@tmpcnt<\@authcnt
   \global\advance\@tmpcnt1
    \raggedright \csname @auth\romannumeral\the\@tmpcnt\endcsname\hfill\\%
   \let\next\@affiliation \vskip1pt
  \else
   \let\next\relax
  \fi
\next}
	 
\newcommand{\affiliation}[1]{%
    \global\advance\@authcnt1
    \expandafter\gdef\csname @auth\romannumeral\the\@authcnt\endcsname
    {#1}}

	
%----------------------------------------------------------------------
%	LIST CONTROL
%----------------------------------------------------------------------  

\RequirePackage{enumitem}
%\setlist{nolistsep} % Uncomment to remove spacing between items in lists (enumerate, itemize)
\setlist{noitemsep}

%----------------------------------------------------------------------
%	ABSTRACT+AUTHOR FRAME
%----------------------------------------------------------------------  

\newcommand{\PaperTitle}[1]{\def\@PaperTitle{#1}}
\newcommand{\Archive}[1]{\def\@Archive{#1}}
\newcommand{\Authors}[1]{\def\@Authors{#1}}
\newcommand{\ExcelYear}[1]{\def\@ExcelYear{#1}}
\newcommand{\Abstract}[1]{\def\@Abstract{#1}}
\newcommand{\Keywords}[1]{\def\@Keywords{#1}}

%   TEASER
\newif\ifHaveTeaser
\HaveTeaserfalse

\newcommand{\Teaser}[1]{\global\HaveTeasertrue\def\@Teaser{#1}}
\newcommand{\TeaserImage}[1]{\includegraphics[height=1in]{#1}}

%   SUPPLEMENTARY
\newif\ifHaveSupplementary
\HaveSupplementaryfalse

\newcommand{\Supplementary}[1]{\global\HaveSupplementarytrue\def\@Supplementary{#1}}


% ---------------------------------------------------------------------

\renewcommand{\@maketitle}{%
\twocolumn[{%
\thispagestyle{empty}%
%\vskip-22pt% {
{\noindent\begin{minipage}[c]{\textwidth}
  \includegraphics[height=34.5mm]{cfde-reports-logo.png}%
  \hfill%
  \raisebox{14mm}{\url{https://cfde.cloud/}}%
\end{minipage}\par}%  
%\includegraphics[width=67mm]{VUT-FIT-logo.pdf}\\[2pt]%
%\includegraphics[width=67mm]{ExcelAtFIT-logo.pdf}%
%\includegraphics[height=14.5mm]{maayanlab-logo.png}\hfill%
%\includegraphics[width=67mm]{ExcelAtFIT-logo.pdf}\\[2pt]%
%\includegraphics[width=67mm]{VUT-FIT-logo.pdf}%
%~~\rotatebox{90}{\color{color-sections}\@ExcelYear}%
\vskip20pt%
{\raggedright\color{color-sections}\sffamily\bfseries\fontsize{20}{25}\selectfont \@PaperTitle\par}%
\vskip10pt%
{\raggedright\color{color-sections}\sffamily\fontsize{12}{16}\selectfont  \@Authors\par}%
\vskip18pt%
\ifHaveTeaser\parbox{\linewidth}{\centering\@Teaser}\vskip10pt\fi%
\fcolorbox{color-blue}{white}{%
\parbox{\textwidth-2\fboxsep-2\fboxrule}{\centering%
\colorbox{color2!8}{%
\parbox{\textwidth-4\fboxsep-2\fboxrule}{%
\sffamily\textbf{\abstractname}\\\@Abstract%\\[4pt]%
%\textbf{\keywordname:} \@Keywords%
%\\[4pt]\textbf{\supplementaryname:} \ifHaveSupplementary\@Supplementary\else\textit{N/A}\fi %
}%
}%
\vskip4pt%
\begingroup%
\raggedright\sffamily\small%
\footnotesize\@affiliation\par%
\endgroup%%
}%
}%
\vskip25pt%
}]%
}%

%%% Make figure placement a little more predictable.
% We trust the user to move figures if this results
% in ugliness.
% Minimize bad page breaks at figures
\renewcommand{\textfraction}{0.01}
\renewcommand{\floatpagefraction}{0.99}
\renewcommand{\topfraction}{0.99}
\renewcommand{\bottomfraction}{0.99}
\renewcommand{\dblfloatpagefraction}{0.99}
\renewcommand{\dbltopfraction}{0.99}
\setcounter{totalnumber}{99}
\setcounter{topnumber}{99}
\setcounter{bottomnumber}{99}




%----------------------------------------------------------------------
%	ADDED STUFF FOR Excel@FIT
%----------------------------------------------------------------------  
\setlength{\columnsep}{0.55cm} % Distance between the two columns of text
\setlength{\fboxrule}{0.75pt} % Width of the border around the abstract

\definecolor{color-sections}{RGB}{3,5,80} % Color of the article title and sections
%\definecolor{color-sections}{RGB}{10,20,120} % Color of the article title and sections
\definecolor{color-links}{RGB}{10,15,160} % Color of the article title and sections
%\definecolor{color-links}{RGB}{170,15,10} % Color of the article title and sections
\definecolor{color-orange}{RGB}{246,135,18} % Color of the article title and sections
\definecolor{color-blue}{RGB}{25,103,179} % Color of the article title and sections
\definecolor{color2}{RGB}{20,20,20} % Color of the boxes behind the abstract and headings
\definecolor{color-ref}{RGB}{25,103,179}
\definecolor{color-lineno}{RGB}{90,140,160}

%\definecolor{color-sections}{RGB}{66,139,202} % Color of the article title and sections


\usepackage[unicode]{hyperref} % Required for hyperlinks
\hypersetup{hidelinks,colorlinks,breaklinks=true,urlcolor=color-links,citecolor=color-ref,linkcolor=color-ref,bookmarksopen=false,pdftitle={Title},pdfauthor={Author},pdfkeywords={Keywords}}

\newcommand{\keywordname}{Keywords} % Defines the keywords heading name
\newcommand{\supplementaryname}{Supplementary Material} % Defines the keywords heading name

\setcounter{tocdepth}{1}

\pagestyle{empty} % Removes page numbering from the first page

\newcommand{\setslovaknames}{%
	\renewcommand{\abstractname}{Abstrakt}%
	\renewcommand{\keywordname}{Kľúčové slová}%
	\renewcommand{\supplementaryname}{Priložené materiály}%
	\renewcommand{\figurename}{Obrázok}%
	\renewcommand{\tablename}{Tabuľka}%
	\renewcommand{\refname}{Literatúra}%
}

\newcommand{\setczechnames}{%
	\renewcommand{\abstractname}{Abstrakt}%
	\renewcommand{\keywordname}{Klíčová slova}%
	\renewcommand{\supplementaryname}{Přiložené materiály}%
	\renewcommand{\figurename}{Obrázek}%
	\renewcommand{\tablename}{Tabulka}%
	\renewcommand{\refname}{Literatura}%
	}

\newcommand{\startdocument}{%
	\ifExcelSlovak\setslovaknames\fi%
	\ifExcelCzech\setczechnames\fi%
	\flushbottom%	
	\maketitle%	
%	\tableofcontents%
	\thispagestyle{empty}%
%	\ifExcelFinal\else\setpagewiselinenumbers\linenumbers\renewcommand\linenumberfont{\scriptsize\color{color-lineno}}\fi%
	}

\renewcommand{\tabcolsep}{1.0mm}
\renewcommand\arraystretch{0.90}

\newif\ifExcelFinal
\ExcelFinalfalse
\def\ExcelFinalCopy{\global\ExcelFinaltrue}

\newif\ifExcelCzech
\ExcelCzechfalse
\def\ExcelCzechLabels{\global\ExcelCzechtrue}

\newif\ifExcelSlovak
\ExcelSlovakfalse
\def\ExcelSlovakLabels{\global\ExcelSlovaktrue}


\@ifclasswith{ExcelAtFIT}{czech}{\global\ExcelCzechtrue}{}
\@ifclasswith{ExcelAtFIT}{slovak}{\global\ExcelSlovaktrue}{}


%--------------------------------------------------------
%	Fix for CZECH and SLOVAK -- needs up to date packages!
%--------------------------------------------------------
\usepackage{regexpatch}
% Change the `-` delimiter to an active character
\@ifclasswith{ExcelAtFIT}{czech}{
	\xpatchparametertext\@@@cmidrule{-}{\cA-}{}{}
	\xpatchparametertext\@cline{-}{\cA-}{}{}
}{}
\@ifclasswith{ExcelAtFIT}{slovak}{
	\xpatchparametertext\@@@cmidrule{-}{\cA-}{}{}
	\xpatchparametertext\@cline{-}{\cA-}{}{}
}{}



%--------------------------------------------------------
%	USEFUL PREPARATIONS
%--------------------------------------------------------

\graphicspath{{images/}{plots/}}

\usepackage{blindtext}
\newcommand{\sgn}{\operatorname{sgn}}

\usepackage{booktabs}

\graphicspath{{images/}{plots/}}

\newcommand{\todo}[1]{\textcolor{red}{\textbf{[[#1]]}}}
\newcommand{\phony}[1]{\textcolor{gray}{#1}}
\newcommand{\blind}[1]{\textcolor{gray}{\Blindtext[#1]}}

%\usepackage[switch]{lineno}

%\setlength{\parindent}{0pt}
\RequirePackage{url}
%\renewcommand{\familydefault}{\sfdefault}

\usepackage[skip=4pt plus3pt, indent=0pt]{parskip}


%----------------------------------------------------------------------
%	REFERENCES
%----------------------------------------------------------------------  

\RequirePackage[
    backend=biber,
    style=nature,
    sorting=none
]{biblatex}
