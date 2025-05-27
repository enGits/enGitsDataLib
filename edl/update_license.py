# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# +                                                                    +
# + This file is part of enGitsDataLib.                                +
# + Copyright 2015-2025 enGits GmbH                                    +
# +                                                                    +
# + enGitsDataLib is released under the MIT License.                   +
# + See LICENSE file for details.                                      +
# +                                                                    +
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

import os
import fnmatch

# List of patterns: [pattern, comment_prefix]
PATTERNS = [
    ('*.cpp', '//'),
    ('*.h', '//'),
    ('*.py', '#'),
    ('CMakeLists.txt', '#'),
]

with open('license_header.txt', 'r', encoding='utf-8') as f:
    LICENSE_HEADER = [line.rstrip() for line in f.readlines()]

for line in LICENSE_HEADER:
    print (line)

def file_matches(filename, pattern):
    return fnmatch.fnmatch(filename, pattern)

def remove_leading_comment_lines(lines, comment_prefix):
    """Remove all consecutive lines at the start that begin with comment_prefix (allowing whitespace before comment)."""
    idx = 0
    prefix_len = len(comment_prefix)
    while idx < len(lines):
        stripped = lines[idx].lstrip()
        if stripped.startswith(comment_prefix):
            idx += 1
        elif stripped == '':  # skip possible blank lines after the header
            idx += 1
        else:
            break
    return lines[idx:]

def insert_license_header(path, comment_prefix):
    with open(path, 'r', encoding='utf-8') as f:
        lines = f.readlines()
    # Remove leading comment lines (and possible blank lines)
    new_content = remove_leading_comment_lines(lines, comment_prefix)
    # Prepare the license block
    license_block = [f"{comment_prefix} {line}\n" for line in LICENSE_HEADER]
    # Insert license header followed by a blank line, then rest of the file
    result = license_block + ['\n'] + new_content
    with open(path, 'w', encoding='utf-8') as f:
        f.writelines(result)
    return True

def main():
    root = os.getcwd()
    updated_files = []
    for dirpath, dirnames, filenames in os.walk(root):
        for pattern, comment_prefix in PATTERNS:
            for filename in filenames:
                if file_matches(filename, pattern):
                    file_path = os.path.join(dirpath, filename)
                    updated = insert_license_header(file_path, comment_prefix)
                    if updated:
                        print(f"Updated: {file_path}")
                        updated_files.append(file_path)
    print(f"Total updated files: {len(updated_files)}")

if __name__ == '__main__':
    main()
