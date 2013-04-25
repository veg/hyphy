from mako.template import Template

ROOT_PATH ="/Users/sweaver/Programming/hyphy/src/"
TAG_FILE = './tags'
TEMPLATE_FILE="./ut_template.cpp"
OUTPUT_DIR = "./"

def get_unique_methods(unique_filename, method_taglist):
    unique_method = filter(lambda x: x[1][0] == unique_filename, method_taglist)
    return list(set([x[0] for x in unique_method]))

def get_includes(fn):
    f = open(fn)
    includes = [x.strip('\n') for x in filter(lambda x: "#include" in x, f)]
    f.close()
    return includes

def get_core_tags_from_file(fn):
    tags_file = open(fn,'r')
    taglist = [tag for tag in tags_file]
    method_taglist = filter(lambda t: "class" in t and "core" in t, taglist)
    method_taglist = [tag.split() for tag in method_taglist]
    method_taglist = [[t[0],(t[1],t[len(t)-1])] for t in method_taglist]
    method_taglist = filter(lambda t: "core" in t[1][0] and ".cpp" in t[1][0], method_taglist)
    tags_file.close()
    return method_taglist

def generate_from_template(fn,class_name,methods,includes):
    new_ut_fn = fn.replace('./core/','./ut_')
    mytemplate = Template(filename=TEMPLATE_FILE)
    msg = mytemplate.render(class_name=class_name,methods=methods,includes=includes)
    f = open(new_ut_fn,'w')
    f.write(msg)
    f.close()


#Get list from file
method_taglist = get_core_tags_from_file(TAG_FILE)

#Get unique filenames
unique_filenames = list(set([t[1] for t in method_taglist]))

#Adapt list to easily view methods by class and file
class_information = [(u[0],u[1].replace("class:",""),
                    get_unique_methods(u[0],method_taglist),
                    get_includes(ROOT_PATH + u[0])) for u in unique_filenames]

#TODO:Don't overwrite existing files
[generate_from_template(info[0],info[1],info[2],info[3]) for info in class_information]
