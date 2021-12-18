import sys
sys.dont_write_bytecode = True 

from helper_sponge_builder import base
from helper_sponge_builder import funcs
import PySimpleGUI as sg
import traceback

sg.theme("SystemDefault")

class GUI():
    pass

class myFrame(sg.Frame):
    def __init__(self, *args, **kw):
        kw["element_justification"] = "left"
        kw["vertical_alignment"] = "top"
        super().__init__(*args, **kw)
        self.sons = {}

class myFrame2(sg.Frame):
    def __init__(self, *args, **kw):
        kw["element_justification"] = "left"
        kw["vertical_alignment"] = "top"
        super().__init__(*args, **kw)
        self.sons = {} 
    def add_row(self,element):
        CurrentRow = self.Rows[0]
        element.Position = (0, len(CurrentRow))
        element.ParentContainer = self
        CurrentRow.append(element)
        if element.Key is not None:
            self.UseDictionary = True
            self.UseDictionary = True

class fatherFrame(sg.Frame):
    def __init__(self, *args, **kw):
        kw["title"] = "Modules"
        kw["vertical_alignment"] = "top"
        kw["layout"]=[[]]
        super().__init__(*args, **kw)
        self.sons = {"核心":myFrame("核心",[[]]), "成键":myFrame("成键",[[]]), "非键":myFrame("非键",[[]]), "环境":myFrame("环境",[[]]), "其他":myFrame("其他",[[]])}
        self.layout([[self.sons["核心"], self.sons["成键"], self.sons["非键"], self.sons["环境"], self.sons["其他"]]])
    def add_row(self,element):
        CurrentRow = self.Rows[1]
        element.Position = (1, len(CurrentRow))
        element.ParentContainer = self
        CurrentRow.append(element)
        if element.Key is not None:
            self.UseDictionary = True


    
def Build_GUI(path):
    gui = GUI()
    gui.Modules = funcs.Search_Modules(path)
    gui.Module_Buttons = []
    gui.Category_Layout = fatherFrame()
    for module in gui.Modules:
        temp = sg.Checkbox(module.Name, tooltip = module.Description, key = module.Name, enable_events = True, 
            default = module.Name in ('common', 'control', 'md_core'), disabled = module.Name in ('common', 'control', 'md_core'))

        gui.Module_Buttons.append(temp)
        temp = module.Category.split("-")
        tempd = gui.Category_Layout
        for i in temp:
            if i not in tempd.sons.keys():
                tempd.sons[i] = myFrame(i,[[]])
                tempd.add_row(tempd.sons[i])
            tempd = tempd.sons[i]
        tempd.add_row(gui.Module_Buttons[-1])
    
    gui.Graph = sg.Multiline("common\ncontrol\nmd_core\n", size = (40, 15), autoscroll=True, key = "GUI_Graph")
    gui.Graph_Contents = []
    gui.Graph_Frame = myFrame("组合框",[[gui.Graph]])
    
    gui.Packages = funcs.Search_Packages(path, gui.Modules)
    gui.Package_Frame = myFrame2("Packages",[[]])
    gui.Package_Buttons = []

    for package in gui.Packages:
        gui.Package_Buttons.append(sg.Button(package.Name, tooltip = package.Description, key = package.Name))
        gui.Package_Frame.add_row(gui.Package_Buttons[-1])
    
    gui.stdout = sg.Multiline(size = (60,10), reroute_stdout = True, autoscroll = True)
    gui.stderr = sg.Multiline(size = (60,10), reroute_stderr = True, autoscroll = True)
    tempButtons = [sg.Button("生成主文件"), sg.Button("生成主头文件"), sg.Button("生成Makefile"), sg.Button("生成图形化界面")]
    if sys.platform == "win32":
        tempButtons.append(sg.Button("直接编译"))
    gui.Operation_Frame = myFrame("操作框",[[myFrame("标准输出",[[gui.stdout]]),myFrame("其他输出",[[gui.stderr]])],tempButtons])
    #gui.Operation_Frame = myFrame("操作框",[[]])
    gui.window = sg.Window("SPONGE",
                [[gui.Category_Layout, sg.VerticalSeparator(), gui.Graph_Frame],[sg.HorizontalSeparator()],[gui.Package_Frame],[sg.HorizontalSeparator()],[gui.Operation_Frame]],
                font = "Arial")
    return gui


if __name__ == "__main__":
    sg.theme("SystemDefault")
    gui = Build_GUI("..")
    while 1:
        event, values = gui.window.read()
        if event == None:
            break
        else:
            try:
                funcs.Event_All(event, values, gui)
            except:
                print(traceback.print_exc(), file=sys.stderr)
    gui.window.close()
