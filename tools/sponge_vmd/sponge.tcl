##
## SPONGE 0.9
##
## A script to visualize SPONGE molecules.
##
## Author: Yijie Xia
##
##  pkg_mkIndex "D:/Program Files/vmd/plugins/noarch/tcl/sponge0.9"

package provide SPONGE 0.9
vmdcon "SPONGE) Loading SPONGE plugin:"


namespace eval SPONGE:: {
    namespace export SPONGE
    
    #创建新分子相关
    variable w_new  ;#创建新文件
    variable mass_in_file_name ;#质量输入文件
    variable bond_in_file_name ;#bond输入文件
    variable crd_in_file_name ; #坐标输入文件
    
    #元素名
    variable name_list {X H He Li Be B C N O F Ne Na Mg Al Si P
        S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As
        Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn
        Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho
        Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po
        At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md
        No Lr Rf Db Sg Bh Hs Mt Ds Rg}
    
    #质量
    variable mass_list {0.00000 1.00794 4.00260 6.941 9.012182 10.811
        12.0107 14.0067 15.9994 18.9984032 20.1797 22.989770
        24.3050 26.981538 28.0855 30.973761 32.065 35.453
        39.948 39.0983 40.078 44.955910 47.867 50.9415
        51.9961 54.938049 55.845 58.9332 58.6934 63.546
        65.409 69.723 72.64 74.92160  78.96 79.904 83.798
        85.4678 87.62 88.90585 91.224 92.90638 95.94 98.0
        101.07 102.90550 106.42 107.8682 112.411 114.818
        118.710 121.760 127.60 126.90447 131.293 132.90545
        137.327 138.9055 140.116 140.90765 144.24 145.0
        150.36 151.964 157.25 158.92534 162.500 164.93032
        167.259 168.93421 173.04 174.967 178.49 180.9479
        183.84 186.207 190.23 192.217 195.078 196.96655
        200.59 204.3833 207.2 208.98038 209.0 210.0 222.0
        223.0 226.0 227.0 232.0381 231.03588 238.02891
        237.0 244.0 243.0 247.0 247.0 251.0 252.0 257.0
        258.0 259.0 262.0 261.0 262.0 266.0 264.0 269.0
        268.0 271.0 272.0}
    
    variable w_add  ;#添加新信息
    variable add_file_name ;#add读入文件
    
    variable backup_frames ;#记录备用的frame
    #备用模块:traj
    namespace eval ReadTraj:: {
        namespace export ReadTraj
        variable box_type ;#记录盒子种类
        variable box_value;#记录盒子信息
        variable box_file ;#记录盒子轨迹文件
    }
}

package provide SPONGE_New 0.9
proc SPONGE::new_button {} {

    variable w_new

    #初始化后就直接打开
    if { [winfo exists .w_new] } {
    wm deiconify $w_new
    return
    }

    #窗口初始化
    set w_new [toplevel ".w_new"]
    wm title $w_new "SPONGE New"
    
    #构局
    #input files frame
    labelframe $w_new.input_frame -text "SPONGE Input Files"
    grid $w_new.input_frame -column 1 -row 1 -columnspan 3
    
    label $w_new.input_frame.first_split -text " "
    grid $w_new.input_frame.first_split 
    
    label $w_new.input_frame.mass_label -text "Mass Input File"
    entry $w_new.input_frame.mass_entry -relief groove -bd 5 -textvariable SPONGE::mass_in_file_name
    button $w_new.input_frame.mass_open -text Open -command SPONGE::New_Button_Mass_Open
    label $w_new.input_frame.mass_hint -text "Required. Mass is used to guess the element."
    label $w_new.input_frame.mass_split -text " "
    
    grid $w_new.input_frame.mass_label -column 1 -row 2 
    grid $w_new.input_frame.mass_entry -column 2 -row 2 
    grid $w_new.input_frame.mass_open -column 3 -row 2 
    grid $w_new.input_frame.mass_hint -column 1 -row 3 -columnspan 3
    grid $w_new.input_frame.mass_split -column 1 -row 4 -columnspan 3
    
    label $w_new.input_frame.bond_label -text "Bond Input File"
    entry $w_new.input_frame.bond_entry -relief groove  -bd 5 -textvariable SPONGE::bond_in_file_name
    button $w_new.input_frame.bond_open -text Open -command SPONGE::New_Button_Bond_Open
    label $w_new.input_frame.bond_hint -text "Optional. This is used to get the bonds"
    label $w_new.input_frame.bond_split -text " "
   
    grid $w_new.input_frame.bond_label -column 1 -row 5
    grid $w_new.input_frame.bond_entry -column 2 -row 5
    grid $w_new.input_frame.bond_open -column 3 -row 5
    grid $w_new.input_frame.bond_hint -column 1 -row 6 -columnspan 3
    grid $w_new.input_frame.bond_split -column 1 -row 7 -columnspan 3
    
    label $w_new.input_frame.crd_label -text "Coordinate Input File"
    entry $w_new.input_frame.crd_entry -relief groove  -bd 5 -textvariable SPONGE::crd_in_file_name
    button $w_new.input_frame.crd_open -text Open -command SPONGE::New_Button_Crd_Open
    label $w_new.input_frame.crd_hint -text "Optional. This is used to get the coordinates"
    label $w_new.input_frame.crd_split -text " "
   
    grid $w_new.input_frame.crd_label -column 1 -row 8
    grid $w_new.input_frame.crd_entry -column 2 -row 8
    grid $w_new.input_frame.crd_open -column 3 -row 8
    grid $w_new.input_frame.crd_hint -column 1 -row 9 -columnspan 3
    grid $w_new.input_frame.crd_split -column 1 -row 10 -columnspan 3
    
    #command frame
    frame $w_new.command_frame
    grid $w_new.command_frame -column 1 -row 2
    
    button $w_new.command_frame.create -text Create -command SPONGE::New_Button_Create
    pack $w_new.command_frame.create -fill x -padx 100

}

proc SPONGE::New_Button_Mass_Open {} {
    variable w_new
    variable mass_in_file_name
    
    set file_types {
    { "Mass Input Files" .txt }
    { "All Files" * }
    }
    
    set mass_in_file_name [tk_getOpenFile -filetypes $file_types \
        -initialdir pwd -initialfile "$mass_in_file_name" \
        -defaultextension .txt]
}

proc SPONGE::New_Button_Bond_Open {} {
    variable w_new
    variable bond_in_file_name
    
    set file_types {
    { "Bond Input Files" .txt }
    { "All Files" * }
    }
    
    set bond_in_file_name [tk_getOpenFile -filetypes $file_types \
        -initialdir pwd -initialfile "$bond_in_file_name" \
        -defaultextension .txt]
}

proc SPONGE::New_Button_Crd_Open {} {
    variable w_new
    variable crd_in_file_name
    
    set file_types {
    { "Coordinate Input Files" .txt }
    { "All Files" * }
    }
    
    set crd_in_file_name [tk_getOpenFile -filetypes $file_types \
        -initialdir pwd -initialfile "$crd_in_file_name" \
        -defaultextension .txt]
}

#参考了topohelpers.tcl by Axel Kohlmeyer <akohlmey@gmail.com>
proc SPONGE::Guess_Name_From_Mass {mass} {
    variable mass_list
    variable name_list
    
    set index 0 ;#找不到设为0
    foreach m $mass_list {
        #直接在列表里找
        if {[expr abs($mass-$m)] < 0.65} {
            set index [lsearch $mass_list $m]
        }
        #氢和氘
        if {($mass > 0.0 && $mass < 2.2)} {
            set index 1
        }
        #铋和钋
        if {($mass > 207.85 && $mass < 208.99)} {
            set index 83
        }
        #镍和钴
        if {($mass > 56.50 && $mass < 58.8133)} {
            set index 27
        }
    }
    set name [lindex $name_list $index]
    return $name
}

#从质量文件创建新分子
proc SPONGE::New_Mol_By_Mass_In_File { mass_in_name } {

    #读质量文件
    set rc [ catch { set fd [open $mass_in_name "r"] } ]
    if { $rc == 1} {
        vmdcon -err "failed to open file $mass_in_name."
        return -1
    }
    #读原子个数
    if { [ gets $fd atom_numbers ] < 0 } {
        vmdcon -err "failed to read atom_numbers in $mass_in_name."
        return -1
    }
    #读原子
    set masses [read $fd]
    set masses [string trim $masses]
    regsub -all {[[:blank:]\n]+} $masses " " masses
    close $fd
    
    
    vmdcon "SPONGE) Atom numbers = $atom_numbers"
    #创建空分子
    set mol -1
    if {[catch {mol new atoms $atom_numbers} mol]} {
        vmdcon -err "failed to create empty molecule: $mol"
        return -1
    }
    set atoms [atomselect $mol all]
    
    #猜测元素符号
    set atom_names {}
    for {set i 0} {$i < $atom_numbers} {incr i} {  
        lappend atom_names [Guess_Name_From_Mass [lindex $masses $i]]
    }
    $atoms set element $atom_names
    $atoms set type $atom_names
    $atoms set name $atom_names
    $atoms set mass $masses
    $atoms delete
    return $mol
}

proc SPONGE::Add_Bond_To_Mol {mol bond_in_name} {
    #读bond信息
    if { [string length $bond_in_name] == 0} {
        return 1
    }
    set rc [ catch { set fd [open $bond_in_name "r"] } ]
    if { $rc == 1} {
        vmdcon -err "failed to open file $bond_in_name."
        return 1
    }
    #键总数与原本的键
    gets $fd bond_numbers
    vmdcon "SPONGE) Bond numbers = $bond_numbers"
    set atoms [atomselect $mol all]
    set atom_numbers [$atoms num]
    for {set i 0} {$i < $atom_numbers} {incr i} {
        set bonds($i) {}
    }
    #读取文件中的键
    for {set i 0} {$i < $bond_numbers} {incr i} {
        gets $fd bond
        set bond [string trim $bond]
        regsub -all {[[:blank:]]+} $bond " " bond
        set bond [split $bond]
        lassign $bond a b k r 
        lappend bonds($a) $b
        lappend bonds($b) $a
    }
    #填回去
    set bonds_new {}
    for {set i 0} {$i < $atom_numbers} {incr i} {
        lappend bonds_new $bonds($i)
    }
    $atoms setbonds $bonds_new
    close $fd
    $atoms delete
    return 0
}

proc SPONGE::Add_Coordinate_To_Mol {mol fname} {
    #打开相关信息
    if { [string length $fname] == 0} {
        return 1
    }
    
    set rc [ catch { set fd [open $fname "r"] } ]
    if { $rc == 1} {
        vmdcon -err "failed to open file $fname."
        return 1
    }
    #获取原子信息，并比较原本原子个数，是否一致
    gets $fd atom_numbers 
    set atoms [atomselect $mol all]
    set atom_numbers_molinfo [$atoms num]
    if { $atom_numbers != $atom_numbers_molinfo } {
        vmdcon -err "The atom_numers in molecule and in coordinate file do not match each other"
        close $fd
        return 1
    }
    #读取坐标信息
    set data {}
    for {set i 0} {$i < $atom_numbers} {incr i} {
        gets $fd crd
        set crd [string trim $crd]
        regsub -all {[[:blank:]]+} $crd " " crd
        set crd [split $crd]
        lappend data $crd
    }
    animate dup $mol
    $atoms set {x y z} $data
    #读取盒子信息
    gets $fd box
    if { [string length $box] == 0} {
        vmdcon "SPONGE) Warning: no box information found in $fname."
    } else {
        pbc set $box
        pbc box
    }
    
    close $fd
    $atoms delete
    
    return 0
}

proc SPONGE::New_Button_Create {} {
    variable w_new
    variable mass_in_file_name
    variable bond_in_file_name
    variable crd_in_file_name
    puts ""
    #读质量文件
    set mol [New_Mol_By_Mass_In_File $mass_in_file_name]
    if {$mol == -1} {
        return 
    }
    #读键连信息
    set rc [Add_Bond_To_Mol $mol $bond_in_file_name]
    if {$rc == 1} {
        vmdcon "SPONGE) Bond Information Not Read"
    }
    
    #读坐标信息
    set rc [Add_Coordinate_To_Mol $mol $crd_in_file_name]
    if {$rc == 1} {
        vmdcon "SPONGE) Coordinate Information Not Read"
    }
    
    mol color [mol default color]
    mol rep [mol default style]
    mol selection [mol default selection]
    mol material [mol default material]
    mol addrep $mol
    display resetview
    mol selupdate 0 $mol on
    return 0
}


proc SPONGE_tk_new {} {
    SPONGE::new_button
    return $SPONGE::w_new
}
vmd_install_extension SPONGE_New SPONGE_tk_new "SPONGE/New Molecule"

package provide SPONGE_Add 0.9

proc SPONGE::add_button {} {

    variable w_add

    #初始化后就直接打开
    if { [winfo exists .w_add] } {
    wm deiconify $w_add
    return
    }

    set w_add [toplevel ".w_add"]
    wm title $w_add "SPONGE Add"
    
    #文件类型下拉框
    frame $w_add.ftype -pady 8
    grid $w_add.ftype -column 1 -row 1
    
    set ftype "trajectory data file"
    label $w_add.ftype.label -text "File Type"
    ttk::combobox  $w_add.ftype.combo -textvariable file_type -stat readonly \
    -values {"trajectory data file" "bond input file" "coordinate input file"}
    
    bind $w_add.ftype.combo <<ComboboxSelected>> { SPONGE::Add_Combobox_Function $file_type}
    
    grid $w_add.ftype.label  -column 1 -row 1
    grid $w_add.ftype.combo  -column 1 -row 2
    
    #主要读取文件
    frame $w_add.fname -pady 8  
    grid $w_add.fname -column 1 -row 2
    #-command SPONGE::New_Button_Add_Open
    label $w_add.fname.label -text "File Name"
    entry $w_add.fname.add_entry -relief groove -bd 5 -textvariable SPONGE::add_file_name
    button $w_add.fname.add_open -text Open -command SPONGE::New_Button_Add_Open
    
    grid $w_add.fname.label -column 1 -row 1 -columnspan 2
    grid $w_add.fname.add_entry -column 1 -row 2 
    grid $w_add.fname.add_open  -column 2 -row 2
    
    #操作按钮
    frame $w_add.command_frame
    grid $w_add.command_frame -column 1 -row 4
    button $w_add.command_frame.load -text Load -command {SPONGE::Add_Button_Load $file_type}
    pack $w_add.command_frame.load -fill x -padx 120
    
    #备用信息
    labelframe $w_add.backup -pady 8 -text "Detailed Options" -labelanchor "n"
    #添加frame后一定要在下面注册一下
    lappend SPONGE::backup_frames "w_add.backup.traj"
    
    
    #读轨迹时候的备用信息
    frame $w_add.backup.traj  
    radiobutton $w_add.backup.traj.radio1 -text "No Box" -variable SPONGE::ReadTraj::box_type -value 1 \
    -command {$SPONGE::w_add.backup.traj.entry2 configure -state disabled; 
              $SPONGE::w_add.backup.traj.entry configure -state disabled; 
              $SPONGE::w_add.backup.traj.open configure -state disabled;}
    radiobutton $w_add.backup.traj.radio2 -text "Fixed Box" -variable SPONGE::ReadTraj::box_type -value 2 \
    -command {$SPONGE::w_add.backup.traj.entry2 configure -state normal;
              $SPONGE::w_add.backup.traj.entry configure -state disabled; 
              $SPONGE::w_add.backup.traj.open configure -state disabled;}
    label $w_add.backup.traj.label2 -text "Format: X Y Z or X Y Z A B C"
    entry $w_add.backup.traj.entry2 -relief groove -bd 5 -textvariable SPONGE::ReadTraj::box_value   
    radiobutton $w_add.backup.traj.radio3 -text "Flexible Box" -variable SPONGE::ReadTraj::box_type -value 3 \
    -command {$SPONGE::w_add.backup.traj.entry2 configure -state disabled;
              $SPONGE::w_add.backup.traj.entry configure -state normal; 
              $SPONGE::w_add.backup.traj.open configure -state normal;}
    label $w_add.backup.traj.label -text "Box Trajectory File"  
    entry $w_add.backup.traj.entry -relief groove -bd 5 -textvariable SPONGE::ReadTraj::box_file
    button $w_add.backup.traj.open -text Open -command {set SPONGE::ReadTraj::box_file [tk_getOpenFile \
        -initialdir pwd -initialfile "$SPONGE::ReadTraj::box_file" \
        -defaultextension *]}
    
    #其他的备用信息
    #暂无
}

proc SPONGE::New_Button_Add_Open {} {
    variable w_add
    variable add_file_name
    
    set file_types {
    { "text Files" .txt }
    { "data Files" .dat }
    { "All Files" * }
    }
    
    set add_file_name [tk_getOpenFile -filetypes $file_types \
        -initialdir pwd -initialfile "$add_file_name" \
        -defaultextension *]
}


proc SPONGE::Add_Combobox_Function {filetype} {
    variable w_add
    
    if { [winfo exists  $w_add.backup] } {
        grid forget  $w_add.backup
    }
    
    foreach path $SPONGE::backup_frames {
        if { [winfo exists $$path] } {
            grid forget $$path
        }
        
    }
    
    if { [string match $filetype "trajectory data file"] } {
        grid $w_add.backup -column 1 -row 3
        grid $w_add.backup.traj -column 1 -row 1
        grid $w_add.backup.traj.radio1 -column 1  -row 1 -sticky W
        grid $w_add.backup.traj.radio2 -column 1  -row 2 -sticky W
        grid $w_add.backup.traj.label2 -column 1  -row 3
        grid $w_add.backup.traj.entry2 -column 1  -row 4
        grid $w_add.backup.traj.radio3 -column 1  -row 6 -sticky W
        grid $w_add.backup.traj.label -column 1 -row 7 -columnspan 2
        grid $w_add.backup.traj.entry -column 1 -row 8 
        grid $w_add.backup.traj.open  -column 2 -row 8
        $w_add.backup.traj.radio1 invoke
    }
}


proc SPONGE::Add_Button_Load {filetype} {
    variable w_add
    variable add_file_name
    if { [string match $filetype "trajectory data file"] } {
        set rc [SPONGE::ReadTraj::function $add_file_name]
        if {$rc == 1} {
            vmdcon "SPONGE) Trajectory Not Read"
        }
    } elseif { [string match $filetype "bond input file"] } {
        set rc [Add_Bond_To_Mol top $add_file_name]
        if {$rc == 1} {
            vmdcon "SPONGE) Bond Information Not Read"
        }
    } elseif { [string match $filetype "coordinate input file"] } {
        set rc [Add_Coordinate_To_Mol top $add_file_name]
        if {$rc == 1} {
            vmdcon "SPONGE) Coordinate Information Not Read"
        }
    }
    
}

proc SPONGE::ReadTraj::function {add_file_name} {
    variable box_type
    variable box_value
    variable box_file
   
    
    set atom_numbers [molinfo top get numatoms]    
    
    if {[catch {open $add_file_name rb} fp]} {
        vmdcon -err "failed to open file $add_file_name"
        return 1
    }
    
    if {$box_type == 3} {
        if {[catch {open $box_file r} fb]} {
            vmdcon -err "failed to open file $box_file"
            return 1
        }
    }
    
    set one_frame_size [expr 12 * $atom_numbers]
    set one_frame_number [expr 3 * $atom_numbers]

        
    set sel [atomselect top all]
    set crds_one_frame_binary [read $fp $one_frame_size]
    while { [string length $crds_one_frame_binary] == $one_frame_size} {
        set data [list]
        binary scan $crds_one_frame_binary f$one_frame_number crd_one_frame
        for {set j 0} {$j < $atom_numbers} {incr j} {
            set crd [lrange $crd_one_frame [expr 3*$j] [expr 3*$j+2]]             
            lappend data $crd 
        }
        animate dup top
        $sel set {x y z} $data
        if {$box_type == 2} {
            pbc set $box_value
        } elseif {$box_type == 3} {
            gets $fb box_value
            pbc set $box_value
        }
        set crds_one_frame_binary [read $fp $one_frame_size]
    }
    display resetview
    if {$box_type > 1} {
        pbc box
    }
    close $fp
    if {$box_type == 3} {
        close $fb
    }
    return 0
}

proc SPONGE_tk_add {} {
    SPONGE::add_button
    return $SPONGE::w_add
}
vmd_install_extension SPONGE_Add SPONGE_tk_add "SPONGE/Add Information"


vmdcon "SPONGE) Finishing SPONGE plugin for VMD"
vmdcon "SPONGE) SPONGE plugin version:"