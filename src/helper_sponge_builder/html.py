from collections import OrderedDict
UNNAMED = -1
class HTML_ELEMENT:
    def __init__(self, html_type, tab, subtab):
        self.type = html_type
        self.contents = None
        self.name = None
        self.nickname = None
        self.tab = tab
        self.subtab = subtab
        self.enable = None
        self.disable = None
        self.default = None
        self.priority = 0
    def __repr__(self):
        return """    name = {}
    nickname = {}
    type = {}
    contents = {}
    tab = {}
    subtab = {}
    enable = {}
    disable = {}
        """.format(*list(map(str,[self.name, self.nickname, self.type, self.contents, self.tab, self.subtab, self.enable, self.disable])))

class HTML_INFO:
    def __init__(self, smf_contents):
        self.elements = []
        current_tab = None
        current_subtab = None
        smf_contents = smf_contents.split(";")
        for smf_content in smf_contents:
            smf_content = smf_content.split("=")
            flag = smf_content[0].strip()
            if flag == "select":
                self.elements.append(HTML_ELEMENT(flag, current_tab, current_subtab))
                self.elements[-1].contents = OrderedDict()
                smf_content = smf_content[1].strip()             
                smf_content = smf_content.split("{")[1:]
                for c in smf_content:
                    c = c.split(":")
                    key = c[0].strip()
                    value = c[1].split("}")[0].strip()
                    self.elements[-1].contents[key] = value
            elif flag == "add_choice":
                self.elements.append(HTML_ELEMENT(flag, current_tab, current_subtab))
                self.elements[-1].contents = smf_content[1].strip()
            elif flag == "input":
                self.elements.append(HTML_ELEMENT(flag, current_tab, current_subtab))
                self.elements[-1].contents = smf_content[1].strip()
            elif flag in ("hint", "up_hint", "down_hint"):
                self.elements.append(HTML_ELEMENT(flag, current_tab, current_subtab))
                self.elements[-1].contents = smf_content[1].strip().replace("\eq","=").replace("\semi",";").replace("\\n", "<br/>").replace("\\t", "&nbsp;"*8).replace(" ", "&nbsp;").replace("\t","&nbsp;"*8).replace("\n", "<br/>")
            elif flag == "name":     
                setattr(self.elements[-1], flag, smf_content[1].strip())
            elif flag == "nickname":
                setattr(self.elements[-1], flag, smf_content[1].strip())
            elif flag == "tab":
                current_tab =  smf_content[1].strip()
                setattr(self.elements[-1], flag, current_tab)
            elif flag == "subtab":
                current_subtab =  smf_content[1].strip()
                setattr(self.elements[-1], flag, current_subtab)
            elif flag == "enable":
                smf_content = smf_content[1].strip()
                self.elements[-1].enable = [smf_content.split()[0].strip()]
                smf_content = smf_content.split("{")[1:]
                value = []
                for c in smf_content:
                    value.append(c.split("}")[0].strip())
                self.elements[-1].enable.append(value)
            elif flag == "disable":
                smf_content = smf_content[1].strip()
                self.elements[-1].disable = [smf_content.split()[0].strip()]
                smf_content = smf_content.split("{")[1:]
                value = []
                for c in smf_content:
                    value.append(c.split("}")[0].strip())
                self.elements[-1].disable.append(value)
            elif flag == "default":
                smf_content = smf_content[1].strip()
                setattr(self.elements[-1], flag, smf_content)
            elif flag == "priority":
                smf_content = float(smf_content[1].strip())
                setattr(self.elements[-1], flag, smf_content)
            elif flag == "":
                continue
            else:
                raise Exception("Unknown flag", flag, smf_contents, "Do you forget to add ';' in the end?")
                
        #检查处理有问题的
        for element in self.elements:
            if not element.name:
                global UNNAMED
                UNNAMED += 1
                element.name = "UNNAME_" + str(UNNAMED)
            if not element.nickname:
                element.nickname = element.name
            if not (element.tab and element.subtab or element.type == "add_choice"):
                raise Exception(element, "has no tab or subtab")


def Html_Script(modules):
    towrite = """
        <script type="text/javascript">
            enable_maps = new Array()
            disable_maps = new Array()\n"""
    
    for module in modules:
        if hasattr(module, "HTML"):
            temp_HTML = module.HTML
        else:
            continue
        SPACES = " " * 12
        for element in temp_HTML.elements:
            if element.enable:
                towrite += """{0}enable_maps.push(new Array("{1}","{2}",{3}))\n""".format(SPACES, element.name, element.enable[0], ",".join(list(map(lambda x:'"'+x+'"', element.enable[1]))))
            if element.disable:
                towrite += """{0}disable_maps.push(new Array("{1}","{2}",{3}))\n""".format(SPACES, element.name, element.disable[0], ",".join(list(map(lambda x:'"'+x+'"', element.disable[1]))))
            
    towrite += "        </script>"
    return towrite

def Html_Tab_Subtab_Contents(modules):
    #写tab部分
    elements = OrderedDict()
    addChoices = OrderedDict()
    for module in modules:
        if hasattr(module, "HTML"):
            temp_HTML = module.HTML
        else:
            continue
        for element in temp_HTML.elements:
            if element.tab not in elements.keys():
                elements[element.tab] = OrderedDict()
            if element.subtab not in elements[element.tab].keys():
                elements[element.tab][element.subtab] = OrderedDict()
            if element.type == "add_choice":
                if element.contents not in addChoices.keys():
                    addChoices[element.contents] = OrderedDict()
                addChoices[element.contents][element.nickname] = element.name
            else:
                elements[element.tab][element.subtab][element.name] = element
    
    for tab in elements.keys():
        for subtab in elements[tab].keys():
            for name, element in elements[tab][subtab].items():
                if element.type == "select" and name in addChoices.keys():
                    element.contents.update(addChoices[element.name])
    del addChoices  
        
    SPACES = " " * 20
    towrite = '                <div class="TAB">'
    temp = '\n                <div class="SUB_TAB">'
    TOTAL_TAB = 0
    TOTAL_SUBTAB = 0
    for i in elements.keys():
        TOTAL_TAB += 1
        towrite += """\n{0}<input id="TAB_{2}" name="TAB" type="radio" class="TAB" onclick="tabshow()"/>
{0}<label for="TAB_{2}" class="TAB">{1}</label>\n""".format(SPACES, i, TOTAL_TAB)
        temp += """\n{0}<div class="TAB_DIV {1}">""".format(SPACES, "TAB_"+str(TOTAL_TAB))
        for j in elements[i].keys():
            TOTAL_SUBTAB += 1
            temp += """\n{0}<input id="SUB_TAB_{3}" name="{2}" type="radio" class="SUB_TAB {2}" onclick="subtabshow()"/>
{0}<label for="SUB_TAB_{3}" class="SUB_TAB">{1}</label>\n""".format(SPACES+"    ", j, "TAB_"+str(TOTAL_TAB), TOTAL_SUBTAB)
        temp += """{0}</div>\n""".format(SPACES)
                
    towrite += '                </div>'
    temp    += '                </div>'
    towrite += temp
    
    TOTAL_TAB = 0
    TOTAL_SUBTAB = 0    
    for tab in elements.keys():
        TOTAL_TAB += 1
        for subtab in elements[tab].keys():
            TOTAL_SUBTAB += 1
            towrite += '\n                <div class="SUB_TAB_DIV TAB_{0} SUB_TAB_{1}">'.format(TOTAL_TAB, TOTAL_SUBTAB)
            for name, element in elements[tab][subtab].items():
                priority = "PRIORITY=%s "%str(element.priority)
                if element.type in ("hint", "up_hint", "down_hint"):
                    if element.default:
                        default = 'for="' + element.default + '" '
                        tooltip = " TOOLTIP"
                    else:
                        default = ""
                        tooltip= ""                    
                    towrite += '\n                    <p id="{0}" {3}{5}class="{1}{4}">{2}</p>\n'.format(name, element.type.upper(), element.contents, default, tooltip, priority)
                elif element.type == "input":
                    if element.default:
                        default = 'value="' + element.default + '" '
                    else:
                        default = ""
                    if element.contents:
                        placeholder = 'placeholder="请输入'+ element.contents + '" '
                    else:
                        placeholder = ""
                    towrite += """\n                    <div {4}class="INPUTS">
                        <label for={0}>{1}</label>
                        <input class="INPUTS" id={0} {2}{3}type="text">
                    </div>\n""".format(name, element.nickname, default, placeholder, priority)
                elif element.type == "select":
                    towrite+= """\n                    <div {2}class="INPUTS">
                        <label for="{0}">{1}</label>
                        <select id="{0}" class="INPUTS">""".format(name, element.nickname, priority)
                    for nickname, name in element.contents.items():
                        if element.default and element.default == nickname:
                            default = 'selected="selected" '
                        else:
                            default = " "
                        towrite += """\n                        <option {2}value ="{0}">{1}</option>""".format(name, nickname, default)
                    towrite+="""\n                        </select>
                    </div>\n"""
            towrite += '                </div>\n'
    
    #写元素部分
    
    return towrite

def Generate_Html(modules):
    towrite = """<!DOCTYPE html>
<html>
    <head>
        <meta charset="UTF-8">
        <meta name="author" content="夏义杰" />
        <title>SPONGE参数生成器</title>
        <link rel="stylesheet" href="GUI_helper/SPONGE.css" media="screen" type="text/css" />
        <script type="text/javascript" src="GUI_helper/jquery.min.js"></script>
        <script type="text/javascript" src="GUI_helper/clipboard.min.js"></script>
        <script type="text/javascript" src="GUI_helper/SPONGE.js"></script>"""

    towrite += Html_Script(modules);

    towrite +="""
    </head>
    <body>
        <div class="BACK_BOARD">
            <h1 class="TITLE">  
                <span>S</span><span>P</span><span>O</span><span>N</span><span>G</span><span>E</span>
            </h1>
            <div class="FRONT_BOARD">
"""

    towrite += Html_Tab_Subtab_Contents(modules);
    
    towrite += """                
            </div>
            <a class="BUTTON" id="NEW_IN_FILE">生成in文件</a>
        </div>
        <div class="OUTPUT">
            <div>
            <a class="BUTTON" id="COPY_BUTTON">复制</a><a class="BUTTON" id="CLOSE_NEW_IN_FILE">关闭</a>
            <br/>
            <pre></pre>
            </div>
        </div>
    </body>
</html>

"""
    return towrite


