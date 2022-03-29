let SessionLoad = 1
let s:so_save = &g:so | let s:siso_save = &g:siso | setg so=0 siso=0 | setl so=-1 siso=-1
let v:this_session=expand("<sfile>:p")
let Tlist_Exit_OnlyWindow =  0 
let Tlist_Show_One_File =  0 
let WebDevIconsNerdTreeGitPluginForceVAlign =  1 
let UltiSnipsRemoveSelectModeMappings =  1 
let WebDevIconsUnicodeDecorateFolderNodesDefaultSymbol = ""
let Tlist_Display_Prototype =  0 
let UltiSnipsEnableSnipMate =  1 
let Tlist_GainFocus_On_ToggleOpen =  0 
let Tlist_Auto_Update =  1 
let Lf_StlColorscheme = "one"
let UltiSnipsDebugHost = "localhost"
let Tlist_File_Fold_Auto_Close =  0 
let UltiSnipsEditSplit = "vertical"
let VM_Cursor_hl = "Cursor"
let WebDevIconsTabAirLineBeforeGlyphPadding = " "
let Tlist_Sort_Type = "order"
let Tlist_Enable_Fold_Column =  1 
let Tlist_Compact_Format =  0 
let UltiSnipsJumpForwardTrigger = "<C-j>"
let VM_Extend_hl = "Visual"
let Tlist_Ctags_Cmd = "ctags"
let DevIconsEnableDistro =  1 
let Taboo_tabs = ""
let WebDevIconsUnicodeDecorateFileNodes =  1 
let UltiSnipsPMDebugBlocking =  0 
let Tlist_Display_Tag_Scope =  1 
let Tlist_Highlight_Tag_On_BufEnter =  1 
let Tlist_Max_Submenu_Items =  20 
let Tlist_WinWidth =  30 
let Tlist_Close_On_Select =  0 
let NERDTreeUpdateOnCursorHold =  1 
let DevIconsDefaultFolderOpenSymbol = ""
let UltiSnipsJumpBackwardTrigger = "<C-k>"
let UltiSnipsDebugPort =  8080 
let Tlist_Show_Menu =  0 
let Tlist_Use_SingleClick =  0 
let UltiSnipsUsePythonVersion =  3 
let UltiSnipsExpandTrigger = "<C-j>"
let DevIconsEnableFoldersOpenClose =  0 
let DevIconsEnableFolderPatternMatching =  1 
let VM_Insert_hl = "Cursor"
let TagList_title = "__Tag_List__"
let WebDevIconsTabAirLineAfterGlyphPadding = ""
let WebDevIconsUnicodeDecorateFileNodesDefaultSymbol = ""
let NERDTreeGitStatusUpdateOnCursorHold =  1 
let WebDevIconsUnicodeDecorateFolderNodesExactMatches =  1 
let WebDevIconsUnicodeByteOrderMarkerDefaultSymbol = ""
let DevIconsAppendArtifactFix =  0 
let WebDevIconsUnicodeDecorateFolderNodes =  1 
let UltiSnipsDebugServerEnable =  0 
let WebDevIconsUnicodeGlyphDoubleWidth =  1 
let Tlist_Max_Tag_Length =  10 
let UltiSnipsListSnippets = "<c-tab>"
let Tlist_Auto_Open =  0 
let DevIconsEnableFolderExtensionPatternMatching =  0 
let Tlist_Auto_Highlight_Tag =  1 
let Tlist_Use_Horiz_Window =  0 
let WebDevIconsNerdTreeAfterGlyphPadding = " "
let Tlist_Process_File_Always =  0 
let ScreenImpl = "Tmux"
let VM_Mono_hl = "Cursor"
let DevIconsArtifactFixChar = " "
let Tlist_WinHeight =  10 
let Tlist_Inc_Winwidth =  1 
let WebDevIconsNerdTreeBeforeGlyphPadding = " "
let WebDevIconsUnicodeDecorateFolderNodesSymlinkSymbol = ""
let Tlist_Use_Right_Window =  0 
silent only
silent tabonly
cd ~/Code/projects/goal-code/src
if expand('%') == '' && !&modified && line('$') <= 1 && getline(1) == ''
  let s:wipebuf = bufnr('%')
endif
set shortmess=aoO
argglobal
%argdel
edit ~/Code/projects/goal-code/scripts/decode/ThetaSeqPlots.jl
argglobal
balt raw.jl
let s:l = 1 - ((0 * winheight(0) + 20) / 41)
if s:l < 1 | let s:l = 1 | endif
keepjumps exe s:l
normal! zt
keepjumps 1
normal! 0
tabnext 1
badd +311 raw.jl
badd +77 utils.jl
badd +172 ~/Code/projects/goal-code/scripts/decode/ThetaSeqPlots.jl
badd +13 ~/Code/projects/goal-code/src/workspace.jl
if exists('s:wipebuf') && len(win_findbuf(s:wipebuf)) == 0 && getbufvar(s:wipebuf, '&buftype') isnot# 'terminal'
  silent exe 'bwipe ' . s:wipebuf
endif
unlet! s:wipebuf
set winheight=1 winwidth=20 shortmess=filnxtToOFAIc
let s:sx = expand("<sfile>:p:r")."x.vim"
if filereadable(s:sx)
  exe "source " . fnameescape(s:sx)
endif
let &g:so = s:so_save | let &g:siso = s:siso_save
set hlsearch
let g:this_session = v:this_session
let g:this_obsession = v:this_session
doautoall SessionLoadPost
unlet SessionLoad
" vim: set ft=vim :
