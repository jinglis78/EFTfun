from ConfigParser import ConfigParser,SafeConfigParser
from collections import OrderedDict
import itertools,re
import Utils

class MyConfigParser(SafeConfigParser):

    def getAndRm(self,section,option):
        val=self.get(section,option)
        self.remove_option(section,option)
        return val

    def getbooleanAndRm(self,section,option):
        val=self.getboolean(section,option)
        self.remove_option(section,option)
        return val

    def itemsAndRm(self,section):
        items=self.items(section)
        self.remove_section(section)
        return items

def expandCfg(cfg_name,overrides):
    replacements={}
    for ovr in overrides:        
        section=ovr.split('=')[0].split('.')[0]
        option=ovr.split('=')[0].split('.')[1]
        value=ovr.split('=')[1]
        if not section in replacements:
            replacements[section]={}
        replacements[section][option]=value
    parser=ConfigParser()
    parser.optionxform = str
    parser.read(cfg_name)
    wildcards={}
    if parser.has_section('Variables'):
        for key,values in parser.items('Variables'):
            wildcards[key]=Utils.split(values)
    sections_new=OrderedDict()
    for section in parser.sections():
        section_new=OrderedDict()
        for option,value in parser.items(section):
            section_new[option]=value
        sections_new[section]=section_new
    for s in replacements:
        for o in replacements[s]:
            print 'will try to override option',o,'in section',s,'in',cfg_name,'with',replacements[s][o]
            if s in sections_new:
                if o in sections_new[s]:
                    print 'replacing old value',sections_new[s][o],'with',replacements[s][o]
                else:
                    print 'adding new option',o,'with value',replacements[s][o]
                sections_new[s][o]=replacements[s][o]
            else:
                print 'adding new section',s,'and new option',o,'with value',replacements[s][o]
                sections_new[s]=OrderedDict()
                sections_new[s][o]=replacements[s][o]

    for section in sections_new:
        for option in sections_new[section]:
            for wc in wildcards:
                if '$('+wc+')' in sections_new[section][option]:
                    value=sections_new[section][option]
                    replacement = ' '.join(wildcards[wc])
                    new_value=value.replace('$('+wc+')',replacement)
                    sections_new[section][option]=new_value
                if '$('+wc+')' in option:
                    replacement = ' '.join(wildcards[wc])
                    new_option=option.replace('$('+wc+')',replacement)
                    if not new_option in sections_new[section]:
                        sections_new[section][new_option]=value



    for section in sections_new:
        for option in sections_new[section]:
            for wc in wildcards:
                if '@('+wc+')' in option:
                    for replacement in wildcards[wc]:
                        new_option=option.replace('@('+wc+')',replacement)
                        if not new_option in sections_new[section]:
                            sections_new[section][new_option]=sections_new[section][option].replace('@('+wc+')',replacement)

    for section in sections_new:
        for option,value in sections_new[section].items():
            for wc in wildcards:
                res=re.findall('@\('+wc+'\d+\)',option)
                if len(res)>0:
                    prod=itertools.product(wildcards[wc],repeat=len(res) )
                    for p in prod:
                        new_option=option
                        new_value=value
                        for i,r in enumerate(res):
                            new_option=new_option.replace(r,p[i])
                            new_value=new_value.replace(r,p[i])
                            if not new_option in sections_new[section]:
                                sections_new[section][new_option]=new_value

    for section in sections_new:
        for option,value in sections_new[section].items():
            for wc in wildcards:
                res=re.findall('@\('+wc+'\d*\)',option)
                if '$('+wc+')' in value or '$('+wc+')' in option or len(res)>0:
                    if option in sections_new[section]:
                        sections_new[section].pop(option)



    config_new=ConfigParser()
    config_new.optionxform = str
    for section in sections_new:
        config_new.add_section(section)
        for key in sections_new[section]:
            config_new.set(section,key,sections_new[section][key])

    config_new.remove_section('Variables')

    new_cfg_name=cfg_name+'.exp'
    with open(new_cfg_name, 'wb') as configfile:
        config_new.write(configfile)
    return new_cfg_name
