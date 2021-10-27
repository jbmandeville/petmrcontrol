#ifndef IMAGEIO
#define IMAGEIO

#include <QMainWindow>
#include <QtConcurrent>
#include <QProgressBar>
#include <QMap>
#include <QString>
#include <QVector>
#include <QVector2D>
#include <QVector3D>
#include <QVector4D>
#include <QFile>
#include <QWidget>
#include <QRunnable>
#include <QDebug>
#include "io.h"

using namespace utilIO;

/*
inline QDataStream &operator>>(QDataStream &in, QVector<short> &vector)
{
   in.readRawData(reinterpret_cast<char*>(&vector), sizeof(vector));
   return in;
}
*/

enum downSampleTypes
{
    downSample1,
    downSample2,
    downSample4,
    downSample21
};
enum CostFunctionChoice
{
    CostFunction_MI,
    CostFunction_WLS,
    CostFunction_NMI,
    CostFunction_NC
};
enum MCAdjustSet
{
    MCAdjustSetRigidBody,
    MCAdjustSetShiftsOnly
};
enum vectorScalars
{
    vectorScalar_Magnitude,
    vectorScalar_x,
    vectorScalar_y,
    vectorScalar_z
};

#define NCOSTBINS 100 // 1% per bin

class bitMap;

class ImageIO
{
    struct NiftiDataType
    {
        int code, size, swapsize;
        QString name;
    };
private:
        int _integerID=-1;     // ID files by an increasing number, in case names are the same
        bVector _loaded;       // [nThreads]; load/reslice/map file in separate threads using initLoaded; setLoaded; isLoaded;
        fileCategory _category=category_none; // categorize datasets for presentation grouping, thresholds, etc.
        int _swap;             // True if byte order on disk is different to CPU
        int _iTime=0;          // points to current volume index
        dPoint3D _xdOrigin={0.,0.,0.};
        iPoint3D _xdDirection={1,1,1};
        vectorScalars _vectorScalar=vectorScalar_Magnitude;

        int extension_type( const QString &extension );
        int read_xdisplay_header();
        int decode_xdisplay_header_keyword( const QString line);
        int readNiftiHeader( bool verbose );

        dMatrix nifti_quatern_to_mat44( float qb, float qc, float qd,
                                        float qx, float qy, float qz,
                                        float dx, float dy, float dz, float qfac );
        float Gauss(float x, float y, float fwhm);
        float inverse_Gauss(float x, float y, float fwhm);
        float delta_filter( float rX );
        float Lanczos_filter( float fx, int iX );
        void nifti_mat44_to_quatern( dMatrix R, float *qb, float *qc, float *qd,
                                     float *qx, float *qy, float *qz,
                                     float *dx, float *dy, float *dz, float *qfac );
        Mat33 nifti_mat33_polar( Mat33 A );
        Mat33 nifti_mat33_inverse( Mat33 R );
        float nifti_mat33_determ( Mat33 R );
        float nifti_mat33_rownorm( Mat33 A );
        float nifti_mat33_colnorm( Mat33 A );

        void swap_16( void *ptr );
        void swap_32( void *ptr );
        dMatrix general_affine_to_mat44( float *srow_x, float *srow_y, float *srow_z );

		static void SwapBytes(size_t n, int siz, void *ar);
		static void SwapNiftiHeader(struct nifti_1_header *h);

        void dilateOrErode(bool dilate, int iTime);

        double power_int(float x, int iPoly);
        double spatialPoly1D(int iPower, int iPos, iPoint2D range);
        double spatialPolyDesign( iPoint3D Dim, iPoint3D iPolyPower, iPoint3D Voxel);

protected:  // can be accessed from derived classes
        QString _imageFileName, _imageHeaderName; // Paths to header and image files
        QString _imageFileIDName;  // a stripped-down name for display: no directory and no extension
        QString _referencedName;  // how the file is referenced at UI (e.g., could include environment variables)
        QString _absoluteFileName; // converted from potentially relative path to absolute path
        void CheckVoxelIndices(fVoxel voxel, const QString callingFunction);
        void CheckVoxelIndices(iPoint3D voxel, const QString callingFunction);
        void CheckVoxelIndices(int iX, int iY, int iZ, const QString callingFunction);
        void CheckVoxelIndices( int i3d, const QString callingFunction );

public:
        imageData _data; // make this public until I figure out how to get it for swapping
        // ImageIO( const ImageIO &imageStack);
        // Public I/O functions
        int readFileHeader(QString fileName, bool verbose );
        int readImageData();
        int readImageData(QProgressBar *progressBar);
        void resliceImageData(ImageIO *templateImage);
        void resliceImageData(ImageIO *templateImage, QProgressBar *progressBar);
        dPoint3D convertVoxelToSpace(iPoint3D iVoxel);
        dPoint3D convertVoxelToSpace( fVoxel voxel);
        dPoint3D convertVoxelToSpace(dPoint3D rVoxel);
        iPoint3D convertSpaceToVoxel(dPoint3D space);
        dPoint3D convertSpaceToVoxelFraction(dPoint3D space);
        double interpolateValueFromSpaceCoord(dPoint3D space, dPoint3D fwhm, iPoint3D wrap,
                                              iPoint3D kernelType, int iTime, bool &inVolume);
        void flattenByGLM(iPoint3D dimPoly, bitMap mask);
        void smoothVolume(int iTime, double fwhm);
        void smoothVolume(int iTime, dPoint3D fwhm);
        void smoothVolumes(double fwhm);
        void smoothVolumes(dPoint3D fwhm);
        void smoothSlices(int iTime, double fwhm);
        void smoothSlices(int iTime, dPoint2D fwhm);
        void smoothSlices(double fwhm);
        void smoothSlices(dPoint2D fwhm);
        void logTransform(bool forward, bitMap mask);
        void logTransform(bool forward);
        inline void dilate(int iTime) {dilateOrErode(true,iTime);}
        inline void erode(int iTime)  {dilateOrErode(false,iTime);}
        inline void close(int iTime)  {dilate(iTime); erode(iTime);}
        inline void open(int iTime)   {erode(iTime); dilate(iTime);}
        double calculateCostFunction(CostFunctionChoice costFunction, bool includeZeroes,
                                     int iVolume, ImageIO *templateVolume, int iVolumeTemplate);
        double getMinDistanceToBoundaryWeight(iPoint3D voxel);
        double calculateAverage(int iTime, ImageIO *templateFile, int iTimeTemplate, bool includeZeroes);
        dPoint2D getExtrema(int iTime);
        dPoint2D getExtrema();
        void resample( iPoint3D dimensionNew, dPoint3D resolutionNew );
        void upSampleOrDownSample(iPoint3D dimensionNew);
        void edgeFilter(int iTime);
        void edgeFilter3D(int iTime);
        void cropUsingOverlay(bitMap mask);

        int writeNiftiFile(QString dirName);
        int writeNiftiData(QString fileName);
        int writeNiftiHeader(QString fileName);

        void setParallelCoordinates(dPoint3D resolution, dPoint3D origin, iPoint3D direction);
        void setParallelCoordinates(dPoint4D resolution, dPoint3D origin, iPoint3D direction);
        void setIjkToXyz(dMatrix ijkToxyz);

        // Public setters
        inline void setCategory(fileCategory category) {_category=category;}
        inline void setVectorScalarMode(vectorScalars scalarMode) {_vectorScalar = scalarMode;}
        inline void setIntegerID(int id) {_integerID = id;}
        inline void setReferencedName(QString fileName) {_referencedName = fileName;}
        inline void setAbsoluteName(QString fileName) {_absoluteFileName = fileName;}
        inline void defineNewMap(QString mapName, fileCategory category, ImageIO *templateImage, int nTime, double DOF)
        { setFileName(mapName); _category=category; copyTemplateDimensions(templateImage, nTime); setDOF(DOF);}
        inline void setMotionCorrectionParameters(int iTime, dVector parameters) {_data.timePoints[iTime].mcParameters = parameters;}

        void setResolution(dPoint4D resolution);
        void setResolution(dPoint3D resolution);
        void setTimeResolution(double dt) {_data.hdr.resolution.t = dt;}
        inline void setStorageFormat(fileStorageTypes iType) {_data.hdr.storageType = iType;}

        void initLoaded (int nThreads);
        void setLoaded();
        void setCurrentVolume(int iTime);
        void setStackValue(const int i3d, const int iT, const double value);
        void setStackValue(const int iX, const int iY, const int iZ, const int iT, const double value);
        void setStackValue(const iPoint3D voxel, const int iT, const double value);
        void setStackValue(const fVoxel voxel, const int iT);
        void setStackValue(const iPoint4D voxel, const double value);

        void setStackValue(const int i3d, const int iT, const Point3D value);
        void setStackValue(const iPoint3D voxel, const int iT, const Point3D value);

        void setParallelCenteredCoordinates();
        void setParallelCenteredCoordinates(iPoint4D dimension, dPoint4D resolution);
        void setParallelCenteredCoordinates(iPoint4D dimension, dPoint3D resolution);
        void setParallelCenteredCoordinates(iPoint3D dimension, dPoint3D resolution);

        void copyTemplate(ImageIO *templateImage);
        void copyTemplate(ImageIO *templateImage, QString name);
        void copyTemplateVolume(ImageIO *templateImage, int iVolume, QString name);
        void copyTemplateVolume(ImageIO *templateImage, int iVolume);
        void copyTemplateDimensions(ImageIO *templateImage, int nTime);
        void setDimensions(const iPoint4D dim);
        void setDimensions(const iPoint3D dim, const int nt);
        void setDimensions(const iPoint3D dim);
        void setDimensions(const int nx, const int nY, const int nZ, const int nT);
        void setOrigin(dPoint3D origin);

        void appendImageFileToTimeData(ImageIO imageFile);
        void deleteFirstVolume();

        inline void setImageHeader(imageHeader *hdr) {_data.hdr = *hdr;}
        void setFileName(QString name);
        void setFileName(QString name, int fileType);
        void setFileIDName(QString ID) {_imageFileIDName = ID;}
        inline void setFilePath(QString fileName) {_absoluteFileName = fileName;}
        inline void setDOF(double dof) {_data.hdr.DOF = dof;}

        // Public Getters
        inline bool isAllocated() {return _data.timePoints.size() != 0;}

        inline QString getFileName() {return _referencedName;}
        inline int getIntegerID() {return _integerID;}
        inline vectorScalars getVectorScalarMode() {return _vectorScalar;}
        inline QString getFileIDName() {return _imageFileIDName;}
        inline int getCurrentVolumeIndex() {return _iTime;}
        inline imageHeader *getImageDataHeader() {return &_data.hdr;}
        inline double getDOF() {return _data.hdr.DOF;}
        bool isLoaded();
        bool noneLoaded();
        inline dVector getMotionCorrectionParameters(int iTime) {return _data.timePoints[iTime].mcParameters;}
        inline iPoint3D getImageDimensions() const
        {
            iPoint3D dim;
            dim.x = _data.hdr.dim.x;  dim.y = _data.hdr.dim.y;  dim.z = _data.hdr.dim.z;
            return dim;
        }
        inline dPoint3D getImageResolution() const
        {
            dPoint3D res;
            res.x = _data.hdr.resolution.x;  res.y = _data.hdr.resolution.y;  res.z = _data.hdr.resolution.z;
            return res;
        }
        inline iPoint4D getDimensions() const { return _data.hdr.dim;}
        inline dMatrix getIjkToXyz() const {return _data.hdr.ijk_to_xyz;}
        inline dMatrix getXyzToIjk() const {return _data.hdr.xyz_to_ijk;}

        inline int nx() const { return _data.hdr.dim.x; }
        inline int ny() const { return _data.hdr.dim.y; }
        inline int nz() const { return _data.hdr.dim.z; }
        inline int nt() const { return _data.timePoints.size(); }
        inline double dx() const { return _data.hdr.resolution.x; }
        inline double dy() const { return _data.hdr.resolution.y; }
        inline double dz() const { return _data.hdr.resolution.z; }
        inline double dt() const { return _data.hdr.resolution.t; }

        inline int index3dFromVoxel(const iPoint3D voxel)
        {
            return voxel.z * _data.hdr.points.y + voxel.y * _data.hdr.points.x + voxel.x;
        }
        inline int index3dFromVoxel(const fVoxel voxel)
        {
            return voxel.z * _data.hdr.points.y + voxel.y * _data.hdr.points.x + voxel.x;
        }
        inline int index3dFromVoxel(const int iX, const int iY, const int iZ)
        {
            return iZ * _data.hdr.points.y + iY * _data.hdr.points.x + iX;
        }

        inline dPoint4D getResolution() const {return _data.hdr.resolution;}
        dPoint3D getOrigin();
        dPoint3D getFOV();
        iPoint3D getDirection();
        iPoint3D determineDirectionMap(dMatrix ijk_to_xyz);
        inline bool hasOffDiagonalTransformation() {return _data.hdr.offDiagonal;}

        inline int voxelsPerSlice() const  { return _data.hdr.points.y; }
        inline int voxelsPerVolume() const { return _data.hdr.points.z; }
        inline int voxelsTotal() const     { return _data.hdr.points.t; }

        inline int datatype() const { return _data.hdr.dataType; }
        inline fileCategory getCategory() {return _category;}
        inline QString getReferencedName() {return _referencedName;}
        inline QString getAbsoluteName() {return _absoluteFileName;}

        dPoint3D coordinateFromVoxel(iPoint3D voxel);
        dPoint3D coordinateFromVoxel(int ix, int iY, int iZ);
        dcVector convertToComplexVolume(int iTime);

        // 4D data getters
        double getStackValue(int iX, int iY, int iZ, int it);
        double getStackValue(int i2d, int iZ, int it);
        double getStackValue(iPoint3D, int it);
        double getStackValue(fVoxel voxel, int it);
        double getStackValue(iPoint4D dim4d);
        double getStackValue(int i3d, int iTime);

        // 3D data getters for volumes with a single time point
        double getVolumeValue(int iX, int iY, int iZ);
        double getVolumeValue(int i3d);
        double getVolumeValue(iPoint3D voxel);
        double getVolumeValue(fVoxel voxel);
        Point3D getVolumeVectorValue(int i3d);
        Point3D getVolumeVectorValue(iPoint3D voxel);

        dPoint3D computeCOM(int iVolume, bool squared, bool inverse);
        dPoint3D computeCOMStdev(int iVolume, bool inverse, dPoint3D COM);
        dPoint2D getMinAndMax(int iTime);
        void fillHistogram(int iTime, double lowerLimit, double upperLimit, bitMap overlay, dVector &histogram);
        double getAverageSignalAfterThreshold(int iTime, double peakFraction);
        void multiplyConstant(int iTime, double multiplier);
};

class bitMap : public ImageIO
{
private:
    VoxelSet _voxelSet;             // (this is a list of voxels, not a volume)
    bool _colorInVoxelSet[NCOLORS]; // keep track of which colors are used to speed display of color overlay maps
    inline bool voxelSurrounded(fVoxel voxel) {iPoint3D voxel3d={voxel.x,voxel.y,voxel.z}; return voxelSurrounded(voxel3d);}
    bool voxelSurrounded(iPoint3D voxel);

public:
    void defineBitMap(iPoint3D dim, dPoint3D res);
    void destroyBitMap();
    void clearBitMap();
    void copyBitMap(bitMap source);

    int readOverlayFile(QString fileName, overlayColor color);
    int readBitmapFile(QString fileName, overlayColor color );
    int writeOverlayFile(QString fileName);
    int writeBitmapFile(QString fileName);
    void setOrAppendBitMap( bitMap overlay, bool clear);
    void setOverlayData(OVLFile ovl, bool clear);
    void setRectangularOverlay( iPoint3D center, iPoint3D width, bool clear);
    void mirrorOverlay();
    void growOverlayField(const iPoint3D seedVoxel, ImageIO *colorFile, double threshold);
    void clusterOverlayField(const iPoint3D seedVoxel);
    iPoint3D getNearestNonzeroVoxel(const iPoint3D seedVoxel);
    void setBitMapVoxelIfTouchingConnectedNeighbor(const fVoxel seedVoxel);
    void dilateOverlayUsingThreshold(ImageIO *file, double thresholdRatio,
                                     bool aboveThreshold, bool twoDimensional, bool fillSurroundedVoxels);
    void fillSurroundedVoxels();
    void dilateOverlay(bool twoDimensional);
    void erodeOverlay(bool twoDimensional);
    inline void dilateOverlay() {dilateOverlay(false);}
    inline void erodeOverlay()  {erodeOverlay(false);}
    inline void openOverlay()  {FUNC_ENTER; erodeOverlay(false); dilateOverlay(false);}
    inline void closeOverlay() {dilateOverlay(false); erodeOverlay(false);}
    inline void openOverlay(bool twoDimensional)  {erodeOverlay(twoDimensional); dilateOverlay(twoDimensional);}
    inline void closeOverlay(bool twoDimensional) {dilateOverlay(twoDimensional); erodeOverlay(twoDimensional);}

    void outerShell(bitMap source);

    void selectOverlayFromThresholdGrayScale(ImageIO *signalFile, bool above, double value);
    void selectOverlayFromThresholdColorScale(ImageIO *signalFile, bool above, double value);
    void applyThresholds(bool reset, ImageIO *sourceFile, int iTime, double lowThresh, double highThresh);

    // Setters
    void setOverlayDimensions(const iPoint3D voxel);
    void setOverlayDimensions(const int nx, const int ny, const int nz);
    void setBitMapValue(const int iX, const int iY, const int iZ, const float value, overlayColor color, bool updateVoxelSet);
    void setBitMapValue(const iPoint3D voxel, const float value, overlayColor color, bool updateVoxelSet);
    void setBitMapValue(const iPoint3D voxel, const float value, bool updateVoxelSet);
    void setBitMapValue(fVoxel ovl, bool updateVoxelSet);
    void setBitMapValueBroadBrush(fVoxel voxel, bool updateVoxelSet);
    void setBitMapValueBroadBrush(const iPoint3D voxel3d, const float value, overlayColor color, bool updateVoxelSet);
    void updateBitMapVoxelSet(overlayColor inputColor, iPoint2D rangeX, iPoint2D rangeY, iPoint2D rangeZ);
    void updateBitMapVoxelSet(overlayColor inputColor);
    inline void updateBitMapVoxelSet() {updateBitMapVoxelSet(Color_Undefined);}
    inline void updateBitMapVoxelSet(iPoint2D rangeX, iPoint2D rangeY, iPoint2D rangeZ)
    { updateBitMapVoxelSet(Color_Undefined, rangeX, rangeY, rangeZ);}
    inline void setFileName(QString fileName) {_imageFileName = fileName; _referencedName = fileName;}
    inline void clearFileName() {_imageFileName = "none"; _referencedName = "none";}
    inline void setColorUsedInOverlay(int iColor, bool state) {_colorInVoxelSet[iColor] = state;}

    // Getters
    inline iPoint3D getDimensions() const {return getImageDimensions();}
    inline dPoint3D getResolution() const {return getImageResolution();}
    inline bool someBitsAreSet() {return (_voxelSet.size() != 0);}
    inline int getBitsSetInBitMap() {return _voxelSet.size();}
    inline VoxelSet getVoxelSet() {return _voxelSet;}
    inline fVoxel getBitMapVoxel(int iVoxel) {return _voxelSet[iVoxel];}
    inline iPoint3D getBitMapVoxel3D(int iVoxel) {fVoxel voxel=_voxelSet[iVoxel]; return {voxel.x,voxel.y,voxel.z};}
    double getNumberVoxelsInOverlay();
    double getVolumeInOverlay();
    overlayColor getBitMapColor(int iX, int iY, int iZ);
    overlayColor getBitMapColor(iPoint3D voxel);
    overlayColor getBitMapColor(fVoxel voxel);
    inline bool colorIsUsedInOverlay(int iColor) {return _colorInVoxelSet[iColor];}
    double getBitMapValue(int i3d);
    double getBitMapValue(int iX, int iY, int iZ);
    double getBitMapValue(iPoint3D voxel);
    double getBitMapValue(fVoxel voxel);
    inline bool isVoxelNonZero(int iVoxel)          {return getBitMapValue(iVoxel)!= 0.f;}
    inline bool isVoxelNonZero(iPoint3D voxel)      {return getBitMapValue(voxel) != 0.f;}
    inline bool isVoxelZero(int iVoxel)             {return getBitMapValue(iVoxel)== 0.f;}
    inline bool isVoxelZero(iPoint3D voxel)         {return getBitMapValue(voxel) == 0.f;}
    inline bool isVoxelAbovePointP5(iPoint3D voxel) {return getBitMapValue(voxel) > 0.5;}
    inline bool isVoxelAbovePointP5(int iVoxel)     {return getBitMapValue(iVoxel) > 0.5;}
    inline QString getFileName() {return _imageFileName;}
    double getMedianSignal(ImageIO *file, bool aboveThreshold, double thresholdRatio);
    iPoint3D computeCOM();
    d3Vector computeCOMs();
};

class voxelAverager : public QObject, public QRunnable
{
    Q_OBJECT   // this can create a problem of "missing vtable"; fix by "run qmake" from Qt Creator build menu
public:
    voxelAverager(int nThreads, VoxelSet voxelList, bVector includeVolume, ImageIO *sourceFile);
    void run(); // if public, it can be run directly in the same thread to test thread safety
private:
    int _nThreads;
    bVector _includeVolume;
    int _nAveraged;
    VoxelSet _voxelList;     // list of slices to map for this thread
    ImageIO *_sourceFile;
signals:
    void progressVoxelAverager(int iProgress);
    void finishedVoxelAverager(VoxelSet voxelList);
public slots:
};

class overlayReader : public QObject, public QRunnable
{
    Q_OBJECT
public:
    OVLFile _file;
    overlayReader(OVLFile overlayFile, iPoint3D dimensions, bool currentSpaceShown);
private:
    void run();
    iPoint3D _dimensions={0,0,0};
    bool _currentSpaceShown;
signals:
    void successfullyReadOverlay(OVLFile imageFile);
    void errorReadingOverlay(QString errorText);
public slots:
};

class voxelReader : public QObject, public QRunnable
{
    Q_OBJECT
public:
    ImageIO _imageFile;
    voxelReader(ImageIO imageFile);
//    int readImageData(QString imageFileName, imageHeader header);
private:
    int machine_byte_order();
    void run();
signals:
    void finishedReadingVoxels(ImageIO imageFile);
    void readProgress(int);
public slots:
//    void readImageData();
};

class voxelSlicer : public QObject , public QRunnable
{
    Q_OBJECT   // this can create a problem of "missing vtable"; fix by "run qmake" from Qt Creator build menu
public:
    voxelSlicer(int firstSlice, int lastSlice, ImageIO imageOriginal, ImageIO imageReslice);
    void run();
private:
    int _iFileType;          // 0 = anatomy, 1 = atlas, 2 = color
    int _firstSlice;         // 1st slice to map
    int _lastSlice;          // last slice to map
    ImageIO  _imageOriginal; // original data that will be resliced
    ImageIO _imageReslice;   // the image data to modify; header must be correct already at initialization
signals:
    void finishedSlice(int whichSlice);
    void finished(int firstSlice, int lastSlice, ImageIO imageFile);
public slots:
};

#endif /* _NIFTI_IO_HEADER_ */
