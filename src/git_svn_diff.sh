#!/bin/csh -f

set PATCH_NAME = "svn"
set FORMAT_PATCH = false
set DO_LOG = false

if ($#argv > 0) then
  foreach arg ( $argv )
    if( $arg == "--format-patch" ) then
      set FORMAT_PATCH = true
    else if( $arg == "--log" ) then
      set DO_LOG = true
    else
      echo "Unknown option '$arg'"
      echo ""
      echo "Available options:"
      echo "  --format-patch \t produce one patch file per commit"
      echo "  --log          \t produce patch summary"
      echo ""
      exit
    endif
  end
endif

# retrieve last svn rev and corresponding git commit
set SVN_REV = `git svn info | grep Revision | cut -d ' ' -f 2`
set GIT_REV = `git svn find-rev r$SVN_REV`

# create diff with git
if ("$FORMAT_PATCH" == "true") then
  set OUT = `git format-patch --no-renames --relative --no-prefix -U0 $GIT_REV`
  set PATCHES = `git log --reverse $GIT_REV.. --oneline | cut -d ' ' -f 1 | sed 's|\(^.*$\)|'$PATCH_NAME'-\1.patch|'`
  set j = 1
  while ( $j <= $#OUT )
    mv $OUT[$j] $PATCHES[$j]
    @ j++
  end
else
  git diff --no-renames --relative --no-prefix -U0 $GIT_REV > $PATCH_NAME.patch
  set PATCHES = $PATCH_NAME.patch
endif

foreach PATCH_FILE ( $PATCHES )
  echo $PATCH_FILE

  # format it to match svn patch
  sed -i '/^new file mode/d' $PATCH_FILE
  sed -i '/^deleted file mode/d' $PATCH_FILE
  sed -i 's|^diff --git \(.*\) .*$|Index: \1|g' $PATCH_FILE
  sed -i 's|^index .*\.\..*$|===================================================================|g' $PATCH_FILE
  sed -i 's|^\(@@ .* @@\).*$|\1|g' $PATCH_FILE
  sed -i 's|^\(--- .*\)$|\1\t(revision '$SVN_REV')|g' $PATCH_FILE
  sed -i 's|^\(+++ .*\)$|\1\t(working copy)|g' $PATCH_FILE
  sed -i 'N;s|^--- \(.*\)\t\(.*\)\n+++ /dev/null.*$|--- \1\t\2\n+++ \1\t(nonexistent)|;P;D' $PATCH_FILE
  sed -i 'N;s|^--- /dev/null.*\n+++ \(.*\)\t\(.*\)$|--- \1\t(nonexistent)\n+++ \1\t\2|;P;D' $PATCH_FILE
end

if ( "$DO_LOG" == "true" ) then
  if ( "$FORMAT_PATCH" == "true" ) then
    git log --reverse --name-status --no-renames --relative --no-prefix $GIT_REV.. > $PATCH_NAME.log
  else
    git diff --name-status --no-renames --relative --no-prefix $GIT_REV > $PATCH_NAME.log
  endif
endif
