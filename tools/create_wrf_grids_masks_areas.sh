#/bin/bash
#
set -u
#
# 1st argument: geogrid filename
geogrid=${1:-geo_em.d01.nc}  # default: geo_em.d01.nc
# 2nd argment: wrf nickname in oasis
model=${2:-atmt}             # default: atmt

# check if $geogrid exists...
if [ ! -f $geogrid ]
then
    echo
    echo "ERROR: file $geogrid not found"
    echo "       please provide the name of the geogrid file when calling $0"
    echo "       e.g.: $0 my_geogrid_file.nc"
    echo
    exit 1
else
  echo "Build grids, masks and areas files based on $geogrid"  
fi

# check if nco functions exist
for nconame in ncwa ncap2 ncks ncrename ncatted
do
    ok=$( which $nconame > /dev/null 2>&1 ; echo $? )
    [ $ok -ne 0 ] && echo "ERROR: $nconame not found" && echo && exit 2
done

# dimension names
xdir=x_$model
ydir=y_$model
cdir=crn_$model

# define ncap2 script
cat <<EOF > ncscript$$
*jj = \$south_north_stag.size ;  // size of y-staggered dimension
*ii = \$west_east_stag.size   ;  // size of y-staggered dimension
defdim("$ydir",jj) ;  // define y-dimension
defdim("$xdir",ii) ;  // define x-dimension
defdim("$cdir", 4) ;  // define corners-dimension
// Warning ncap2 does not accept variables with a dot in their name
// grids variables
lon[          \$${ydir},\$${xdir}] = 0.f ;  // default center lon definition
lat[          \$${ydir},\$${xdir}] = 0.f ;  // default center lat definition
clo[\$${cdir},\$${ydir},\$${xdir}] = 0.f ;  // default corner lon definition
cla[\$${cdir},\$${ydir},\$${xdir}] = 0.f ;  // default corner lat definition
lon(  0:jj-2,0:ii-2) = XLONG_M(0,0:jj-2,0:ii-2) ;  // M-center lon
clo(0,0:jj-2,0:ii-2) = XLONG_C(0,0:jj-2,0:ii-2) ;  // bottom-left  corner lon
clo(1,0:jj-2,0:ii-2) = XLONG_C(0,0:jj-2,1:ii-1) ;  // bottom-right corner lon
clo(2,0:jj-2,0:ii-2) = XLONG_C(0,1:jj-1,1:ii-1) ;  // upper-right  corner lon
clo(3,0:jj-2,0:ii-2) = XLONG_C(0,1:jj-1,0:ii-2) ;  // upper-left   corner lon
lat(  0:jj-2,0:ii-2) = XLAT_M( 0,0:jj-2,0:ii-2) ;  // M-center lat
cla(0,0:jj-2,0:ii-2) = XLAT_C( 0,0:jj-2,0:ii-2) ;  // bottom-left  corner lat
cla(1,0:jj-2,0:ii-2) = XLAT_C( 0,0:jj-2,1:ii-1) ;  // bottom-right corner lat
cla(2,0:jj-2,0:ii-2) = XLAT_C( 0,1:jj-1,1:ii-1) ;  // upper-right  corner lat
cla(3,0:jj-2,0:ii-2) = XLAT_C( 0,1:jj-1,0:ii-2) ;  // upper-left   corner lat
// mask variable
msk[\$${ydir},\$${xdir}] = 1 ;                     // default mask definition (i.e. not active)
msk(0:jj-2,0:ii-2) = LANDMASK(0,0:jj-2,0:ii-2) ;   // fill mask values of M-points
// areas variable
srf[\$${ydir},\$${xdir}] = 1.f ;                   // default surface definition
srf(0:jj-2,0:ii-2) = ( @DX * @DY ) / MAPFAC_MX(0,0:jj-2,0:ii-2) / MAPFAC_MY(0,0:jj-2,0:ii-2) ;
EOF
ncap2 -O -v -S ncscript$$ $geogrid tmp_gma$$.nc          # apply ncap2 script

# grids.nc
ncks -O -v lon,lat,clo,cla tmp_gma$$.nc grids.nc         # keep only the grids variables
ncrename -v lon,$model.lon -v lat,$model.lat \
	 -v clo,$model.clo -v cla,$model.cla grids.nc    # rename variables
ncatted -O -h -a ".*,global,d,," grids.nc                # clean all global attributes
echo "  grids.nc done"

# masks.nc
ncks -O -v msk tmp_gma$$.nc masks.nc        # keep only the masks variable
ncrename -v msk,$model.msk masks.nc         # rename variable
ncatted -O -h -a ".*,global,d,," masks.nc   # clean all global attributes
echo "  masks.nc done"

# areas.nc
ncks -O -v srf tmp_gma$$.nc areas.nc        # keep only the areas variable
ncrename -v srf,$model.srf areas.nc         # rename variable
ncatted -O -h -a ".*,global,d,," areas.nc   # clean all global attributes
echo "  areas.nc done"

# cleaning of temporary files
rm -f ncscript$$ tmp_gma$$.nc     
