#! /bin/csh -f
#
# The purpose of this script is to test changes you have made to your
# local SVN repository by shipping those changes to the Uintah
# buildbot server to run tests against.
#
#
# Before you can run this script, the "buildbot" package must be
# installed on your computer.
#
# - On Debian (as root) you can run: apt-get install buildbot
#
#
# Other caveats:
#
# - You can't run the try if you have new files in your tree...  The
#     only way around this is to "svn commit" the new files first.
#
#
# This is an example of a successful submission:
#
#   % ./buildbot_try_test.sh
#   2018-04-24 21:52:50-0600 [-] Log opened.
#   2018-04-24 21:52:50-0600 [-] using 'pb' connect method
#   2018-04-24 21:52:52-0600 [-] job created
#   2018-04-24 21:52:52-0600 [-] Starting factory <twisted.spread.pb.PBClientFactory instance at 0x106f2cf38>
#   2018-04-24 21:52:52-0600 [Broker,client] Delivering job; comment= None
#   2018-04-24 21:52:52-0600 [Broker,client] job has been delivered
#   2018-04-24 21:52:52-0600 [Broker,client] not waiting for builds to finish
#   2018-04-24 21:52:52-0600 [-] Stopping factory <twisted.spread.pb.PBClientFactory instance at 0x106f2cf38>
#   2018-04-24 21:52:52-0600 [-] Main loop terminated.

#
# Buildbot test script. Execute this script in the top level source dir:
# i.e. /Uintah/trunk/src

# 1. Before running this script the src code must be fully up to date.

# 2. If you are adding new files they must be checked in first.

#______________________________________________________________________
#

set trunkServers = ("Trunk:opt-full-try" "Trunk:dbg-full-try" "Trunk:opt-gpu-try")

set BUILDERS = ""
set CREATE_PATCH = false
set MY_PATCH     = false
set PATCH_FILE = "buildbot_patch.txt"
set USE_GIT = `svn info >& /dev/null && echo false || git svn info >& /dev/null && echo true || echo false`

set VC = "--vc=svn"
set REBOSITORY = ""
set DIFF = ""

# No args so all tests

if ($#argv == 0) then
  foreach server ($trunkServers )
    set BUILDERS = "$BUILDERS --builder=$server"
  end
endif

#__________________________________
# parse inputs
while ( $#argv )
  #echo "($1)"

  # remove punctuation chars and convert to lowercase
  set arg = `echo $1 | tr '[:upper:]' '[:lower:]'`

  switch ($arg:q)
    case createpatch:
      set CREATE_PATCH = true
      shift
      breaksw

    case mypatch:
      set MY_PATCH = true
      set PATCHFILE = $2
      shift; shift
      breaksw

    case trunk-opt:
      set BUILDERS = "$BUILDERS --builder=Trunk:opt-full-try"
      shift
      breaksw
    case trunk-d*b*g:         # debug or dbg
      set BUILDERS = "$BUILDERS --builder=Trunk:dbg-full-try"
      shift
      breaksw

    case trunk-gpu:
      set BUILDERS = "$BUILDERS --builder=Trunk:opt-gpu-try"
      shift
      breaksw
    
    case kokkos-opt:
      set BUILDERS = "$BUILDERS --builder=Kokkos:opt-full-try"
      shift
      breaksw
      
    case all:
      foreach server ($trunkServers)
        set BUILDERS = "$BUILDERS --builder=$server"
      end  
      shift
      breaksw

    default:
      echo " Error parsing inputs."
      echo " Usage: buildbot_try.sh   [options]"
      echo "             Options:"
      echo "              trunk-opt                  Trunk:opt-full-try server"
      echo "              trunk-debug/dbg            Trunk:dbg-full-try server"
      echo "              trunk-gpu                  Trunk:opt-gpu-try server"
      echo "              kokkos-opt                 Kokkos:opt-full-try server"
      echo "              all                        run trunk(opt + dbg + gpu) try servers"
      echo "              create_patch               run svn diff on src/ and submit that patch"
      echo "              myPatch      <patchFile>  submit the patchfile to the try servers"
      echo "   Now exiting"
      exit(1)
      breaksw
  endsw
end

#______________________________________________________________________
#
# Note normally svn will automatically create a patch.  It will
# contain context lines which are superfluous and make the patch
# bigger than necessary.

# If your changes are huge manually create an svn patch with no
# context - still has a max 640 Mbytes.

# always create patch
if ( $USE_GIT == "true" ) then
  /bin/rm -rf "$PATCH_FILE" >& /dev/null

  # retrieve last svn rev and corresponding git commit
  set SVN_REV = `git svn info | grep Revision | cut -d ' ' -f 2`
  set GIT_REV = `git svn find-rev r$SVN_REV`
  # create diff with git (ommitting deleted as svn diff)
  git diff --no-renames --relative --no-prefix -U0 $GIT_REV > $PATCH_FILE
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

  ls -l $PATCH_FILE

  set VC = ""
  set DIFF = "--diff=$PATCH_FILE"
  set REPOSITORY = "--repository=https://gforge.sci.utah.edu/svn/uintah/trunk/src"

else if( $CREATE_PATCH == "true" ) then
  /bin/rm -rf "$PATCH_FILE" >& /dev/null

  svn diff -x --context=0 > $PATCH_FILE

  ls -l $PATCH_FILE

  set VC = ""
  set DIFF = "--diff=$PATCH_FILE"
  set REPOSITORY = "--repository=https://gforge.sci.utah.edu/svn/uintah/trunk/src"
endif
#__________________________________
# use a user created patch

if( $MY_PATCH == "true" ) then
  if( ! -e $PATCHFILE ) then
    echo "  Error:  Could not find the patch file $PATCHFILE"
    exit 1
  endif

  set PATCH = "--diff=$PATCHFILE --repository=https://gforge.sci.utah.edu/svn/uintah/trunk/src"
endif

#__________________________________

echo "  PATCH $PATCH"
echo "  BUILDERS: $BUILDERS"

#__________________________________

echo buildbot --verbose try \
         --connect=pb \
         --master=uintah-build.chpc.utah.edu:8031 \
         --username=buildbot_try \
         --passwd=try_buildbot \
         --topdir=. \
         --who=`whoami` \
         $VC $DIFF $REPOSITORY $BUILDERS

#buildbot --verbose try \
#         --connect=pb \
#         --master=uintah-build.chpc.utah.edu:8031 \
#         --username=buildbot_try \
#         --passwd=try_buildbot \
#         --topdir=. \
#         --who=`whoami` \
#         $VC $DIFF $REPOSITORY $BUILDERS

echo $status

# cleanup
if( $CREATE_PATCH == "true" || $USE_GIT == "true" ) then
  /bin/rm -rf $PATCH_FILE >& /dev/null
endif

exit
