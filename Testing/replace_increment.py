import sys
import re

def replace_increments(filepath):
    with open(filepath, 'r') as f:
        content = f.read()
    
    # ONLY + is allowed! mi = mi + 1 -> mi += 1
    pattern = r'\b([a-zA-Z_][a-zA-Z0-9_]*(?:\[[^\]]+\])*)\s*=\s*\1\s*\+\s*([^;]+?)\s*;'
    
    def replacer(match):
        var = match.group(1)
        val = match.group(2)
        return f"{var} += {val};"
    
    new_content = re.sub(pattern, replacer, content)
    
    pattern_for = r'\b([a-zA-Z_][a-zA-Z0-9_]*(?:\[[^\]]+\])*)\s*=\s*\1\s*\+\s*([^;)]+?)\s*(\)|;)'
    def replacer_for(match):
        var = match.group(1)
        val = match.group(2)
        end = match.group(3)
        return f"{var} += {val}{end}"
        
    new_content = re.sub(pattern_for, replacer_for, new_content)
    
    with open(filepath, 'w') as f:
        f.write(new_content)

if __name__ == "__main__":
    replace_increments(sys.argv[1])
