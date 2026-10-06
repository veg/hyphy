import sys

def format_braces(filepath):
    with open(filepath, 'r') as f:
        lines = f.readlines()
        
    out_lines = []
    i = 0
    while i < len(lines):
        line = lines[i]
        stripped = line.strip()
        
        if stripped == "{":
            # If the current line is just "{", we should append it to the previous line
            if out_lines:
                prev_line = out_lines.pop()
                # Remove trailing newline and add " {"
                out_lines.append(prev_line.rstrip() + " {\n")
        elif stripped.startswith("{") and len(stripped) > 1:
            # If the line starts with "{" but has other stuff, e.g. "{ return 1; }"
            # Maybe just append " {" and the rest of the line on a new line, or same line?
            # Usually it's just "{" on a line by itself in HyPhy. Let's just handle "{" by itself for now.
            if out_lines:
                prev_line = out_lines.pop()
                out_lines.append(prev_line.rstrip() + " {\n")
                out_lines.append(line.replace("{", "", 1).lstrip())
            else:
                out_lines.append(line)
        else:
            out_lines.append(line)
        i += 1
        
    with open(filepath, 'w') as f:
        f.writelines(out_lines)

if __name__ == "__main__":
    format_braces(sys.argv[1])
