; Based on the example from http://www.seas.gwu.edu/~drum/java/lectures/appendix/installer/install.html

!define VERSION "2.2.4"
!define PACKAGE_NAME "HyPhy"

Icon "../../src/gui/res/Windows/desk.ico"
UninstallIcon "../../src/gui/res/Windows/desk.ico"

XPStyle on
ShowInstDetails hide
ShowUninstDetails hide

; Name of our application
Name "${PACKAGE_NAME} ${VERSION}"

; The file to write
OutFile "${PACKAGE_NAME}${VERSION}.exe"

; Set the default Installation Directory
InstallDir "$PROGRAMFILES\${PACKAGE_NAME} ${VERSION}"

; Set the text which prompts the user to enter the installation directory
DirText "Please choose a directory to which you'd like to install this ${PACKAGE_NAME} ${VERSION}."

; ----------------------------------------------------------------------------------
; *************************** Pages *******************************
; ----------------------------------------------------------------------------------

Page license
Page directory
Page components
Page instfiles

UninstPage uninstConfirm
UninstPage instfiles

; ----------------------------------------------------------------------------------
; *************************** Lincense *******************************
; ----------------------------------------------------------------------------------

!define licensefile license.rtf
!ifdef licensefile
LicenseText "License"
LicenseData "${licensefile}"
!endif


; ----------------------------------------------------------------------------------
; *************************** SECTION FOR INSTALLING *******************************
; ----------------------------------------------------------------------------------

AutoCloseWindow false
ShowInstDetails show

Section "Create Desktop Shortcut" SectionX
    SetShellVarContext current
    CreateShortCut "$DESKTOP\HyPhy.lnk" "$INSTDIR\HYPHY.exe"
SectionEnd

Section "Install HyPhy" ; A "useful" name is not needed as we are not installing separate components

; Set output path to the installation directory. Also sets the working
; directory for shortcuts

SetOutPath $INSTDIR\

File ../../HYPHY.EXE
File /r ../../help/*.pdf
File /r ../../res/*.*

WriteUninstaller $INSTDIR\Uninstall.exe

; ///////////////// CREATE SHORT CUTS //////////////////////////////////////

CreateDirectory "$SMPROGRAMS\${PACKAGE_NAME} ${VERSION}"
CreateShortCut "$SMPROGRAMS\${PACKAGE_NAME} ${VERSION}\Run ${PACKAGE_NAME} ${VERSION}.lnk" "$INSTDIR\HYPHY.exe"
CreateShortCut "$SMPROGRAMS\${PACKAGE_NAME} ${VERSION}\Getting Started.lnk" "$INSTDIR\Getting Started With HyPhy.pdf"
CreateShortCut "$SMPROGRAMS\${PACKAGE_NAME} ${VERSION}\Selection Analyses.lnk" "$INSTDIR\SelectionAnalyses.pdf"
CreateShortCut "$SMPROGRAMS\${PACKAGE_NAME} ${VERSION}\Uninstall ${PACKAGE_NAME} ${VERSION}.lnk" "$INSTDIR\Uninstall.exe"

; ///////////////// END CREATING SHORTCUTS //////////////////////////////////

; //////// CREATE REGISTRY KEYS FOR ADD/REMOVE PROGRAMS IN CONTROL PANEL /////////

WriteRegStr HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\${PACKAGE_NAME} ${VERSION}" "DisplayName"\
"${PACKAGE_NAME} ${VERSION} (remove only)"

WriteRegStr HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\${PACKAGE_NAME} ${VERSION}" "UninstallString" \
"$INSTDIR\Uninstall.exe"

; //////////////////////// END CREATING REGISTRY KEYS ////////////////////////////

MessageBox MB_OK "Installation was successful. Please use the shortcuts in the 'Start' menu to launch HyPhy and check its documentation."

SectionEnd

; ----------------------------------------------------------------------------------
; ************************** SECTION FOR UNINSTALLING ******************************
; ----------------------------------------------------------------------------------

Section "Uninstall"
; remove all the files and folders
Delete $INSTDIR\Uninstall.exe ; delete self
Delete $INSTDIR\HYPHY.exe
RMDIR /r $INSTDIR

RMDir $INSTDIR

; now remove all the startmenu links
Delete "$SMPROGRAMS\${PACKAGE_NAME} ${VERSION}\Run ${PACKAGE_NAME} ${VERSION}.lnk"
Delete "$SMPROGRAMS\${PACKAGE_NAME} ${VERSION}\Uninstall ${PACKAGE_NAME} ${VERSION}.lnk"
Delete "$SMPROGRAMS\${PACKAGE_NAME} ${VERSION}\Getting Started.lnk"
Delete "$SMPROGRAMS\${PACKAGE_NAME} ${VERSION}\SelectionAnalyses.lnk"
Delete "$DESKTOP\HyPhy.lnk"

RMDIR /r "$SMPROGRAMS\${PACKAGE_NAME} ${VERSION}"

; Now delete registry keys
DeleteRegKey HKEY_LOCAL_MACHINE "SOFTWARE\${PACKAGE_NAME} ${VERSION}"
DeleteRegKey HKEY_LOCAL_MACHINE "SOFTWARE\Microsoft\Windows\CurrentVersion\Uninstall\${PACKAGE_NAME} ${VERSION}"

SectionEnd
