#!/usr/bin/env python3

## \file org2dox.py
# \brief A simple utility for converting README.org into .dox format.
# \author Kevin K. Chen (Princeton University)
#
# See the file's docstring for a complete description of this file.

## \namespace org2dox
# \brief Python namespace for \c org2dox.py.

"""A simple utility for converting README.org into .dox format.

After reading text from the standard input, the following changes are made:

- "=code=" is changed to "`code`".

- "/** \mainpage" is appended to the beginning of the text, and "*/" is appended
  to the end.

- Instances of "* [...]" at the beginning of lines are changed to "\section
  [...] [...]".  Note that if [...] contains spaces, then only the first word is
  preserved in the first output; e.g.,

  * Additional remarks

  becomes

  \section Additional Additional remarks

- The above is also done with "**" and "\subsection", and "***" and
  "\subsubsection".

- Lines with "#+BEGIN_SRC" and "#+END_SRC" are deleted.

The result is then printed to the standard output.

"""

import argparse
import re
import sys


## Parse the command-line arguments, process the text, and print it.
def main():
    parse()

    text = sys.stdin.read()
    text = process(text)
    print(text)


## Parse the input arguments for the help flag.
def parse():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.parse_args()


## Make the substitutions in the text.
# \param[in] text The text to process.
# \return The processed text.
def process(text):
    # Start the comment block and add the mainpage command.
    text = '/**\n\mainpage ' + text
    # Turn =code fragments= into `code fragments`.
    text = re.sub(r"""([\s'"({]|^)=(.*?\S)=([\s,;:.?!'"\\)}]|$)""", r'\1`\2`\3',
                  text)
    # Reformat the sections.
    text = re.sub(r'^\* +(\w+)(.*)', r'\section \1 \1\2', text, flags=re.M)
    text = re.sub(r'^\*\* +(\w+)(.*)', r'\subsection \1 \1\2', text, flags=re.M)
    text = re.sub(r'^\*\*\* +(\w+)(.*)', r'\subsubsection \1 \1\2', text,
                  flags=re.M)
    # Remove Org-mode source block delimiters.
    text = re.sub(r'^\s*#\+(BEGIN|END)_SRC.*$', r'', text, flags=re.M)
    # End the comment block.
    text += '*/'

    return text


if __name__ == '__main__':
    main()
