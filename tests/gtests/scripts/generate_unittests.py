from mako.template import Template
from itertools import chain
from functools import lru_cache
import create_constructors
import re

ROOT_PATH ="/Users/sweaver/Programming/hyphy/src/"
TAG_FILE = './tags'
TEMPLATE_FILE="./gtest.cpp.tpl"
OUTPUT_DIR = "../core/"
CONSTRUCT_DICT = create_constructors.get_constructs()

def get_constructors_from_file(fn):
    tags_file = open(fn,'r')
    taglist = [tag for tag in tags_file]
    method_taglist = filter(lambda t: "class" in t and "core" in t, taglist)
    method_taglist = [tag.split('\t') for tag in method_taglist]
    method_taglist = [[t[0],(get_return_type(t[2]),get_parameters(t[2])),(t[1],t[len(t)-1])] for t in method_taglist]
    method_taglist = list(filter(lambda t: t[0].strip()==t[len(t)-1][1].replace("class:","").strip(), method_taglist))
    tags_file.close()
    # Convert to dictionary
    new_dict = {}
    for constructor in method_taglist:
        new_dict[constructor[0]]=constructor[1:]
    return new_dict

def get_constructor_params(classname):
    constructor_list = get_constructors_from_file(TAG_FILE)
    if classname in constructor_list.keys():
        return constructor_list[classname][0][1]
    else:
        return None

#Get return types of method in ctags
def get_return_type(definition):
    try:
        results = re.search('\^ *([\w\*].*?) (\**)[\w:]+\s*\(',definition)
        rs = ''.join(results.groups()).replace(" ","")
        rs = rs.replace("virtual ","")
        rs = rs.replace("unsigned","")
        rs = rs.replace("const","")
        return rs
    except:
        try:
            results = re.search('\^(.*?) ',definition).group(1)
            return results
        except:
            pass
        return ''

#Get parameters of method in ctags
def get_parameters(definition):
    try:
        param_list = re.search('\((.*?)[\)\?\$]',definition).group(1).split(',')
        param_list = list(filter(lambda x:x!='',param_list))
        return param_list
    except:
        return None

def get_unique_methods(unique_filename, method_taglist):
    unique_method = filter(lambda x: x[2][0] == unique_filename and "::" not in x[1][0], method_taglist)
    unique_list = list(unique_method)
    return unique_list

def get_dependent_classes(name,c):
    #We need to return both the class, and the parameters the class needs
    total_params = list(set(list(chain.from_iterable([x[1][1] if x[1][1]!=None else [] for x in c]))))
    total_params = filter(lambda x: x != "" and x != "...", list(set(total_params)))
    unique_types = []
    for t in total_params:
        split_t = t.split()
        if len(split_t) > 1:
            param = split_t[0].strip()
        else:
            param = t.strip()
        if "FILE" not in param:
            unique_types.append(param)
    unique_types = [t.replace("&","") for t in unique_types]
    unique_types = [t.replace("*","") for t in unique_types]
    unique_types = filter(lambda x: x != "const", list(set(unique_types)))
    unique_types = filter(lambda x: x != "void", list(set(unique_types)))
    tuples = [(u,CONSTRUCT_DICT[u][0][0]) if u in CONSTRUCT_DICT.keys() else (u,None) for u in unique_types]
    try:
        tuples.append((name,CONSTRUCT_DICT[name][0][0]))
        depends = CONSTRUCT_DICT[name][0][1]
        try:
            [tuples.append((name,CONSTRUCT_DICT[name][0][0])) for name in depends]
        except:
            pass
    except:
        tuples.append((name,name + "test = new " + name + "()"))
    return list(set(tuples))

def get_includes(fn):
    f = open(fn)
    includes = [x.strip('\n') for x in filter(lambda x: x.strip().startswith("#include"), f)]
    includes = list(filter(lambda x: "dmalloc" not in x, includes))
    includes = list(filter(lambda x: "timer" not in x, includes))
    includes = list(filter(lambda x: "profiler" not in x, includes))
    includes = list(filter(lambda x: "HYTreePanel.h" not in x, includes))
    includes = list(filter(lambda x: "HYChartWindow.h" not in x, includes))
    includes = list(filter(lambda x: "HYUtils.h" not in x, includes))
    includes = list(filter(lambda x: "Windows.h" not in x, includes))
    includes = list(filter(lambda x: "windows.h" not in x, includes))
    includes = list(filter(lambda x: "HYSharedMain.h" not in x, includes))
    includes = list(filter(lambda x: "HYConsoleWindow.h" not in x, includes))
    includes = list(filter(lambda x: "opencl_kernels.h" not in x, includes))
    includes = list(filter(lambda x: "ocl" not in x, includes))
    includes = list(filter(lambda x: "OpenCL" not in x, includes))
    includes = list(filter(lambda x: "opencl" not in x, includes))
    includes = list(filter(lambda x: "preferences" not in x, includes))
    f.close()
    return includes

def get_core_tags_from_file(fn):
    tags_file = open(fn,'r')
    taglist = [tag for tag in tags_file]
    method_taglist = filter(lambda t: "class" in t and "core" in t, taglist)
    method_taglist = [tag.split('\t') for tag in method_taglist]
    method_taglist = [[t[0],(get_return_type(t[2]),get_parameters(t[2])),(t[1],t[len(t)-1])] for t in method_taglist]
    method_taglist = list(filter(lambda t: "core" in t[2][0] and ".cpp" in t[2][0], method_taglist))
    method_taglist = list(filter(lambda t: "file" not in t[2][1], method_taglist))
    tags_file.close()
    return list(method_taglist)

def generate_from_template(fn,class_name,methods,includes,dependencies):
    new_ut_fn = fn.replace('./core/','./' + OUTPUT_DIR + 'ut_')
    if "ut_strings" in new_ut_fn:
        return
    if "ut_calcnode" in new_ut_fn:
        return
    mytemplate = Template(filename=TEMPLATE_FILE)
    msg = mytemplate.render(class_name=class_name,methods=methods,includes=includes,objects=dependencies)
    f = open(new_ut_fn,'w')
    f.write(msg)
    f.close()

def stringify_method_parameters(method):
    param_list = method
    if param_list == None or len(param_list) == 0 or param_list[0] == '':
        return ''
    type_list = []
    for param in param_list:
        split_types = param.split()
        split_types = list(filter(
                            lambda x: "const" not in x and "unsigned" not in x and "void" not in x,
                            split_types))
        if len(split_types)==0:
            return ''
        type = split_types[0] + "test"
        if len(split_types) > 1:
            if not split_types[1].startswith("*"):
                type = "*" + type
        else:
            type = "*" + type
        type_list.append(type)
    return ", ".join(type_list)

#Get list from file
method_taglist = get_core_tags_from_file(TAG_FILE)

#Get unique filenames
unique_filenames = list(set([t[2] for t in method_taglist]))

#Adapt list to easily view methods by class and file
class_information = [{"filename"     : u[0],
                      "classname"    : u[1].replace("class:","").strip(),
                      "methods"      : get_unique_methods(u[0],method_taglist),
                      "includes"     : get_includes(ROOT_PATH + u[0]),
                      "dependencies" : get_dependent_classes(u[1].replace("class:","").strip(),get_unique_methods(u[0],method_taglist))}
                      for u in unique_filenames]


#Operator list
#arithmetic operators: + - * / % and += -= *= /= %= (all binary infix); + - (unary prefix); ++ -- (unary prefix and postfix)
#bit manipulation: & | ^ << >> and &= |= ^= <<= >>= (all binary infix); ~ (unary prefix)
#boolean algebra: == != < > <= >= || && (all binary infix); ! (unary prefix)
#memory management: new new[] delete delete[]
#implicit conversion operators
#miscellany: = [] -> , (all binary infix); * & (all unary prefix) () (function call, n-ary infix)

operators =      [",","+","-","*","/","%","+=","-=","*=","/=","%=","+","-","++","--","&","|","^","<<",">>","&=","|=","^=","<<=",">>=","~","==","!=","<",">","<=",">=","||","&&","!","->","*","()","[]","="]
operator_names = ["Comma","Plus","Subtract","Mult","Div","Mod","PlusEqual","SubEqual","MultEqual","DivEqual","ModEqual","UnaryPlus","UnaryNeg","Inc","Dec","Amp","Pipe","Carot","DoubleLess","DoubleGreater","AmpEqual","PipeEqual","CarotEqual","DoubleLessEqual","DoubleGreaterEqual","Tilde","DoubleEqual","NotEqual","LessThan","GreaterThan","LessThanEqual","GreaterThanEqual","DoublePipe","DoubleAmp","Not","PointerClass","Pointer","Parenths","Brackets","Equal"]
operator_dict = dict(zip(operators,operator_names))

#We need to adapt the methods to have proper parameters
for c in class_information:
    count = 1
    methods = c["methods"]
    for key in range(len(methods)):
        #if "file" in methods[key][2][1]:
        #    new_tuple = (methods[key][2][0], "_Matrix")
        #    methods[key][2] = new_tuple
        #    c["classname"] = "_Matrix"
        if "~" in methods[key][0]:
            methods[key][0] = methods[key][0].replace("~","Deconstructor")
        if "operator " in methods[key][0]:
            try:
                methods[key][0] = "operator" + operator_dict[methods[key][0].split()[1]]
            except:
                methods[key][0] = "unknownOperator" + "".join(methods[key][0].split())
        if key > 0:
            if methods[key-count][0] == methods[key][0]:
                methods[key].append(methods[key][0] + str(count))
                count+=1
            else:
                methods[key].append(methods[key][0])
                count=1
        else:
            methods[key].append(methods[key][0])

        method_strings = stringify_method_parameters(methods[key][1][1])
        methods[key].append(method_strings)

#TODO:Don't overwrite existing files
[generate_from_template(info["filename"],info["classname"],info["methods"],info["includes"],info["dependencies"]) for info in class_information]
