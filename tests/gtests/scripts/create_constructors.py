import re

TEMPLATE_FILE="./gtest_constructor.cpp.tpl"
TAG_FILE = './tags'

#Get parameters of method in ctags
def get_parameters(definition):
    try:
        return re.search('\((.*?)\)',definition).group(1).split(',')
    except:
        return None

def get_return_type(definition):
    try:
        return re.search('\^(.*?) ',definition).group(2)
    except:
        return None

def clean_params(l):
  param_list = l
  if param_list == None or param_list[0] == '':
      return ''
  type_list = []
  for param in param_list:
      split_types = param.split()
      split_types = list(filter(
                          lambda x: "const" not in x and "unsigned" not in x and "void" not in x,
                          split_types))
      if len(split_types)==0:
          return ''
      type = split_types[0]
      if len(split_types) > 1:
          if not split_types[1].startswith("*"):
              type = type
      else:
          type = type
      type_list.append(type)
      return type_list

def stringify_method_parameters(method):
    param_list = method
    if param_list == None or param_list[0] == '':
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

#We need to create constructors of each kind of object
def get_constructors_from_file(fn):
    tags_file = open(fn,'r')
    taglist = [tag for tag in tags_file]
    method_taglist = filter(lambda t: "class" in t and "core" in t, taglist)
    method_taglist = [tag.split('\t') for tag in method_taglist]
    method_taglist = [[t[0],(get_return_type(t[2]),t[2]),(t[1],t[len(t)-1])] for t in method_taglist]
    method_taglist = filter(lambda t: t[0].strip()==t[len(t)-1][1].replace("class:","").strip(), method_taglist)
    tags_file.close()
    method_taglist = [(t[0],get_parameters(t[1][1])) for t in method_taglist]
    return list(method_taglist)

#Get list from file

def construct(name,params):
  if params == None:
      params = ''
  return name + "test = new " + name + "(" + params + ")"


def get_constructs():
    ct = get_constructors_from_file(TAG_FILE)
    #Create a dictionary with object name and instiate strings
    construct_dictionary = {}
    for t in ct:
        if t[0] not in construct_dictionary.keys():
            construct_dictionary.update({t[0]:[(construct(t[0],stringify_method_parameters(t[1])),clean_params(t[1]))]})
        else:
            construct_dictionary[t[0]].append((construct(t[0],stringify_method_parameters(t[1])),clean_params(t[1])))
    return construct_dictionary

