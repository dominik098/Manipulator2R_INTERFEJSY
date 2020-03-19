%% Modified from olplainarticle.cls
% A plain single-column layout for articles.
%
% Created by Overleaf (writeLaTeX).
%
% Based on the SelfArx document class.
%
\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{physiome}[2020/03/13, v1.0]
\RequirePackage[utf8]{inputenc}
\RequirePackage[T1]{fontenc}
\RequirePackage{textcomp}
\RequirePackage[english]{babel}

\RequirePackage{etoolbox}
\RequirePackage{xparse}
\RequirePackage{calc}
\AtEndOfClass{\RequirePackage{microtype}}
\newif\if@reqslineno
\DeclareOption{lineno}{\@reqslinenotrue}
\DeclareOption*{\PassOptionsToClass{\CurrentOption}{article}}
\ProcessOptions*
\LoadClass{article}
\def\@journallogo{Physiome-under}
\RequirePackage{mathpazo}
\RequirePackage{helvet}
\RequirePackage{lineno}     % For numbering the lines if lineno option is given
\if@reqslineno\linenumbers\fi

\RequirePackage{ifpdf}

\RequirePackage{amsmath,amssymb}
\RequirePackage{graphicx,xcolor}
\RequirePackage{booktabs}
\RequirePackage{authblk}
\RequirePackage{mdframed}
\RequirePackage{marginnote}

\RequirePackage[left=5cm,%
                right=2cm,%
                top=2.25cm,%
                bottom=2.25cm,%
                headheight=12pt,%
                letterpaper,%
                reversemarginpar,%
                marginparwidth=3cm,%
                marginparsep=2em%
                ]{geometry}%
                
\RequirePackage[labelfont={bf,sf},%
                labelsep=period,%
                justification=raggedright]{caption}
                
\RequirePackage[round]{natbib}
\RequirePackage{multibib}
\newcites{primepub}{Primary Publications}
\bibliographystyle{unsrtnat}
\bibliographystyleprimepub{unsrtnat}

\newcommand{\primarypubs}[2]{%
  \nociteprimepub{#2}
    \begin{mdframed}[backgroundcolor=color2!50!white,topline=false,bottomline=false,rightline=false,leftline=false,font=\sffamily,leftmargin=-10\p@,rightmargin=-10\p@,innertopmargin=-5\p@,innerbottommargin=-5\p@]
  \if@reqslineno\internallinenumbers\fi
  \renewcommand{\bibsection}{\subsection*{\refname}}
  \bibliographyprimepub{#1}
  \end{mdframed}\vskip18pt
}

%
% headers and footers
%
\RequirePackage{fancyhdr}  % custom headers/footers
\RequirePackage{lastpage}  % Number of pages in the document
\pagestyle{fancy}          % Enables the custom headers/footers
% Headers
\lhead{}%
\chead{}%
\rhead{}%
% Footers
\lfoot{}%
\cfoot{}%
\rfoot{\small\sffamily\bfseries\thepage/\pageref{LastPage}}%
\renewcommand{\headrulewidth}{0pt}% % No header rule
\renewcommand{\footrulewidth}{0pt}% % No footer rule

%
% section/subsection/paragraph set-up
%
\RequirePackage[explicit]{titlesec}
\titleformat{\section}
  {\color{color1}\large\sffamily\bfseries}
  {\thesection}
  {0.5em}
  {\MakeUppercase{#1}}
  []
\titleformat{name=\section,numberless}
  {\color{color1}\large\sffamily\bfseries}
  {}
  {0em}
  {\MakeUppercase{#1}}
  []  
\titleformat{\subsection}
  {\sffamily\bfseries}
  {\thesubsection}
  {0.5em}
  {#1}
  []
\titleformat{\subsubsection}
  {\sffamily\small\bfseries\itshape}
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
\titlespacing*{\subsection}{0pc}{2.5ex \@plus3pt \@minus2pt}{3pt}
\titlespacing*{\subsubsection}{0pc}{2ex \@plus2.5pt \@minus1.5pt}{2pt}
\titlespacing*{\paragraph}{0pc}{1.5ex \@plus2pt \@minus1pt}{10pt}

%
% tableofcontents set-up
%
\RequirePackage{titletoc}
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
  
\RequirePackage{enumitem}
%\setlist{nolistsep} % Uncomment to remove spacing between items in lists (enumerate, itemize)


%
% article meta data
%
\newcommand{\keywords}[1]{\noindent\textcolor{color1}{\keywordname}\quad #1\par}%\vskip25\p@}

%
% custom title page 
%
\newcommand{\publicationdate}[1]{\def\physiome@pubdate{#1}}
\newcommand{\curator}[1]{\def\physiome@curator{#1}}
\let\authblk@author\author
\RenewDocumentCommand{\author}{ o o m }{%
  \authblk@author[#1]{#3%
    \IfValueT{#2}{\thanks{Corresponding Author: \detokenize{#2}}}%
  }%
}
\renewcommand\AB@authnote[1]{\textsuperscript{#1}}
\renewcommand\AB@affilnote[1]{\textsuperscript{#1}}
\renewcommand{\@maketitle}{%
    \marginnote{\includegraphics[width=\hsize]{\@journallogo}%
    \ifdefvoid{\physiome@pubdate}{}{\\\parbox{\hsize}{\centering\footnotesize Published \physiome@pubdate}}%
    }%
    {\raggedright\sffamily\bfseries\fontsize{20}{25}\selectfont\@title\par}%
    \vskip10pt
    {\raggedright\@author\par}
    \vskip18pt%
}
\apptocmd{\maketitle}{\thispagestyle{empty}}{}{}
%-----------------------------------------------

\renewenvironment{abstract}{%
  \section*{\abstractname}%
  \begin{mdframed}[topline=false,bottomline=false,rightline=false,leftline=false,font=\sffamily,leftmargin=-10\p@,rightmargin=-10\p@]
  \if@reqslineno\internallinenumbers\fi
}{%
  \end{mdframed}\vskip18pt
}


\setlength{\columnsep}{0.55cm} % Distance between the two columns of text
% \setlength{\fboxrule}{0.75pt} % Width of the border around the abstract

\definecolor{color1}{RGB}{0,0,0} % Color of section headings
\definecolor{color2}{RGB}{209,178,220} % Color of the box behind the abstract
\newcommand{\keywordname}{Keywords:} % Defines the keywords heading name

\renewcommand\Authfont{\fontsize{12}{14}\sffamily\bfseries}
\renewcommand\Affilfont{\fontsize{10}{12}\sffamily\bfseries}

\newlength{\tocsep} 
\setlength\tocsep{1.5pc} % Sets the indentation of the sections in the table of contents
\setcounter{tocdepth}{3} % Show only three levels in the table of contents section: sections, subsections and subsubsections

\RequirePackage{lipsum} % Required to insert dummy text
%-----------------------------------------------
\renewcommand{\bibsection}{%
  \ifdefvoid{\physiome@curator}{}{%
  \bigskip
  \noindent\emph{Curated by \physiome@curator}\par
  }
  \section{\refname}
}
