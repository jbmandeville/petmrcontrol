#ifndef IO_H
#define IO_H

#include <QString>
#include <QVector>
#include <QProcess>
#include <QtMath>
#include <QDebug>

#define FUNC_ENTER qDebug() << Q_FUNC_INFO << "enter"
#define FUNC_EXIT  qDebug() << Q_FUNC_INFO << "exit"
#define FUNC_INFO  qDebug() << Q_FUNC_INFO
#define FUNC_ERROR qDebug() << "Error in" << Q_FUNC_INFO

#define PI 3.141592654
#define PI_180 (PI/180.)

#define NCOLORS 14 // do not include undefined value
enum overlayColor
{   //                 R   G   B
    Color_Cyan,    //  0  255 255
    Color_Teal,    //  0  128 128
    Color_Purple,  // 128 128  0
    Color_Magenta, // 255  0  255
    Color_Pink,    // 255 192 203
    Color_Red,     // 255  0   0
    Color_Lime,    //  0  255  0
    Color_Green,   //  0  128  0
    Color_Sky,     //  0  191 255
    Color_Blue,    //  0   0  255
    Color_Yellow,  // 255 255  0
    Color_Orange,  // 255 165  0
    Color_Tan,     // 210 180 140
    Color_Brown,   // 139 69  19
    Color_Undefined // this must be last, so that the enum values match the combobox
};

enum fileStorageTypes
{
    BFLOAT    = 0,
    BSHORT    = 1,
    BLONG     = 2,
    UINT8     = 3,
    USHORT    = 4,
    BDOUBLE   = 5,
    NIFTI_NII = 6,
    NIFTI_IMG = 7
};

enum fileCategory
{
    category_none,
    category_timeSeries,
    category_anatomy,
    category_atlas,
    category_summary,
    category_preProcessingSource,
    category_color,
    category_pValue,
    category_TValue,
    category_signalValue,
    category_BPValue,
    category_Relaxation,
    category_ADC,
    category_anisotropy,
    category_vector,
    category_complex,
    fileCategory_transient // files generated within the program but not saved
};

struct iPoint2D {union {int x; int lower;}; union {int y; int upper;}; };
struct iPoint3D {int x, y, z;};
struct iPoint4D {int x, y, z, t;};
struct bPoint3D {bool x, y, z;};
struct Point2D  {float x, y;};
struct Point3D  {float x, y, z;};
struct Point4D  {float x, y, z, t;};
struct dPoint2D {union {double x; double lower;}; union {double y; double upper;}; };
struct dPoint3D {double x, y, z;};
struct dPoint4D {double x, y, z, t;};
struct Mat33 {double m[3][3];};
struct Mat44 {double m[4][4];};
struct fComplex {float real, imag;};
struct dComplex {double real, imag;};
typedef QVector<int>                iVector;
typedef QVector<bool>               bVector;
typedef QVector<QChar>              cVector;
typedef QVector<QString>            sVector;
typedef QVector<double>             dVector;
typedef QVector<float >             fVector;
typedef QVector<iPoint2D>          i2Vector;
typedef QVector<iPoint3D>          i3Vector;
typedef QVector<dPoint2D>          d2Vector;
typedef QVector<dPoint3D>          d3Vector;
typedef QVector<fComplex>          fcVector;
typedef QVector<dComplex>          dcVector;
typedef QVector<Point3D>            f3Vector;
typedef QVector<QVector<QString>>   sMatrix;
typedef QVector<QVector<bool>>      bMatrix;
typedef QVector<QVector<QChar>>     cMatrix;
typedef QVector<QVector<int>>       iMatrix;
typedef QVector<QVector<float>>     fMatrix;
typedef QVector<QVector<double>>    dMatrix;
typedef QVector<QVector<fComplex>> fcMatrix;
typedef QVector<QVector<Point3D>>  f3Matrix;
typedef QVector<QVector<dPoint2D>> d2Matrix;
typedef QVector<QVector<dPoint3D>> d3Matrix;
typedef QVector<QVector<iPoint2D>> i2Matrix;
typedef QVector<QVector<QVector<int>>>    iMatrix3;
typedef QVector<QVector<QVector<double>>> dMatrix3;
struct fVoxel {int x, y, z; float value; overlayColor color;};  // can hold a voxel location, value, and optionally a color for display
typedef QVector<fVoxel> VoxelSet;

struct OVLFile
{  // typical overlay-list format: [name] [path] [optional color]
    iPoint3D dimensions={0,0,0};  // keeps image dimensions (not always available)
    QString name;                  // e.g. putamen
    QString path;                  // e.g. $TemplateDir/putamen.ovl
    overlayColor color = Color_Cyan;  // Attach a color to an overlay file
    int iAtlasRegionID;  // indicates region from atlas, with ID as given (e.g region 21); set to 0 for non-atlas regions
    VoxelSet voxelList;// vector of (x,y,z,value)
    bool saved=true;
    bool loaded=false;
};

struct imageHeader
{ // all info to interpret a file should be stored in the header in order to facilitate voxel-reading in a different thread
    QString fileName;       // original file name (if read from file), including directory path
    iPoint4D dim={0,0,0,0}; // # voxels in x, y, z, t
    dPoint4D resolution={0.,0.,0.,0.};    // resolutions in x, y, z
    iPoint4D points={0,0,0,0};        // (nx, nx*ny, nx*ny*nz, nx*ny*nz*nt)
    dMatrix ijk_to_xyz;     // 4x4 transformation from data indices (i,j,k) to space (x,y,z)
    dMatrix xyz_to_ijk;     // 4x4 inverse transformation
    int byte_order;         // byte order 0 = Mac, 1 = PC
    bool offDiagonal=false ;// if true, the transformation matrix does not have parallel voxels & coordinates
    int fileType=NIFTI_NII; // e.g., bshort, bfloat, blong, nifti_nii, nifti_img
    int dataType=0;         // 0 = float, 1 = real-imaginary, 2 = magnitude-phase, 3 = vector
    int storageType=BFLOAT; // storage type: bshort, bfloat, blong
    int slice_code = 0;
    float slice_duration=0.;
    double DOF=0.;          // degrees of freedom for statistics; default value 0
    int voxelOffset;        // offset to the first voxel data in the file
    float scaleSlope;       // scaling factor for data (data --> slope*data + offset); ignore if 0; default 0
    float scaleOffset;      // offset for data; default 0
};

struct timePointData
{
    fVector  f1;            // 3-dimensional volume data vector [i3d]
    fcVector fc;            // 3-dimensional complex vector
    f3Vector f3;            // 3-dimensional 3-vector
    dVector mcParameters={0.,0.,0.,0.,0.,0.};   // motion-correction parameters
};

// Define a structure for a stack of images
struct imageData
{
    imageHeader hdr;
    QVector<timePointData> timePoints;  // allocated as a time vector per file
    iVector color;                      // 3-dimensional color (for bitmaps/overlays)
    bool empty=false;                   // if true, all data is zero

//    iPoint4D dim={0,0,0,0}; // # voxels in x, y, z, w=t
//    Point4D resolution;    // resolutions in x, y, z, t
//    Point4D origin;        // origins in x, y, z, t
//    Point3D direction;     // directions in x, y, z
//    iPoint4D points;       // (nx, nx*ny, nx*ny*nz, nx*ny*nz*nt)
//    Mat44 ijk_to_xyz;      // transformation from data indices (i,j,k) to space (x,y,z)
//    Mat44 xyz_to_ijk;      // inverse transformation
};

struct ROI_data
{
    QString fileName;// Vector of file names for each run in the ROI
    QString name;    // some identifier
    iPoint3D voxel;  // voxel for 1-point ROI
    double nspace;   // # points in space (non-integer due to potential non-binary voxel weights)
    double dt;       // time step
    double dof;
    dVector xTime;   // 1-dimensional data vector [itime]
    dVector ySignal; // 1-dimensional data vector [itime]
};

struct MCPackage // motion-correction package
{
    int iFile;                                  // file index
    int iTime;                                  // volume index
    dVector MCParameters={0.,0.,0.,0.,0.,0.};   // x,y,z,phiX,phiY,phiZ
};

namespace utilIO
{
    template <class T> T* allocate_vector(int n)
    {
        T* a = new T[n];
        return a;
    }
    template <class T> void delete_vector(T *vec)
    {
        delete vec;
    }

    int machine_byte_order();
    void swap(short *ptr);
    void swap(unsigned short *ptr);
    void swap(float *ptr);
    void swap(int  *ptr);
    void swap_16( void *ptr );
    void swap_32( void *ptr );

    inline int index3d(iPoint3D dim, iPoint3D voxel) {return voxel.z*dim.x*dim.y + voxel.y*dim.x + voxel.x;}
    inline int index3d(iPoint3D dim, int iX, int iY, int iZ) {return iZ*dim.x*dim.y + iY*dim.x + iX;}

    int readOverlayFile(VoxelSet &voxelList, QString fileName, overlayColor color, iPoint3D dimSpace);

    void delayMS( int millisecondsToWait );

    QString readTimeTableFile(QString fileName, QStringList &columnNames, dMatrix &table);

}

namespace tk
{
// band matrix solver
class band_matrix
{
private:
    dMatrix m_upper;  // upper band
    dMatrix m_lower;  // lower band
public:
    // constructors/destructors
    band_matrix() {};
    inline band_matrix(int dim, int n_u, int n_l) {resize(dim, n_u, n_l);}
    ~band_matrix() {};                       // destructor
    void resize(int dim, int n_u, int n_l);  // init with dim,n_u,n_l
    inline int dim() const {if ( m_upper.size() > 0) return m_upper[0].size(); else return 0;}
    int num_upper() const
    {
        return m_upper.size()-1;
    }
    int num_lower() const
    {
        return m_lower.size()-1;
    }
    // access operator
    double & operator () (int i, int j);            // write
    double   operator () (int i, int j) const;      // read
    // we can store an additional diogonal (in m_lower)
    // second diag (used in LU decomposition), saved in m_lower
    inline double saved_diag(int i) const {Q_ASSERT( (i>=0) && (i<dim()) ); return m_lower[0][i];}
    inline double &saved_diag(int i) {Q_ASSERT( (i>=0) && (i<dim()) ); return m_lower[0][i];}

    void lu_decompose();
    dVector r_solve(const dVector& b) const;
    dVector l_solve(const dVector& b) const;
    dVector lu_solve(const dVector& b,bool is_lu_decomposed=false);
};

// spline interpolation
class spline
{
public:
    enum bd_type {
        first_deriv = 1,
        second_deriv = 2
    };

private:
    dVector m_x,m_y;            // x,y coordinates of points
    // interpolation parameters
    // f(x) = a*(x-x_i)^3 + b*(x-x_i)^2 + c*(x-x_i) + y_i
    dVector m_a,m_b,m_c;        // spline coefficients
    double  m_b0, m_c0;                     // for left extrapol
    bd_type m_left, m_right;
    double  m_left_value, m_right_value;
    bool    m_force_linear_extrapolation;

public:
    // set default boundary condition to be zero curvature at both ends
    spline(): m_left(second_deriv), m_right(second_deriv),
        m_left_value(0.0), m_right_value(0.0),
        m_force_linear_extrapolation(false)
    {
        ;
    }

    // optional, but if called it has to come be before set_points()
    void set_boundary(bd_type left, double left_value,
                      bd_type right, double right_value,
                      bool force_linear_extrapolation=false);
    void set_points(const dVector& x,
                    const dVector& y, bool cubic_spline=true);
    double operator() (double x) const;
};

}

namespace utilMath
{
#define SQR(a) ((a) == 0.0 ? 0.0 : (a)*(a))
#define NMAT 3

    dComplex complexMultiply(dComplex dc1, dComplex dc2, bool star);
    double computePhase(dComplex dc);
    double computeAmplitude(dComplex dc);
    double dotProduct(dVector v1, dVector v2);
    double dotProduct(dPoint3D v1, dPoint3D v2);
    void FFTW1D_Volume( iPoint3D dim, dcVector &volume, int iDimFFT, bool forward );

    Mat44 invertMat44(Mat44 matrix44 );
    bool dInvertSquareMatrix(dMatrix &dSquareMatrix );

    double deltaKernel(double x);
    double LanczosKernel(double x, double width);
    double GaussKernel(double x, double y, double fwhm);
    double inverse_Gauss(double x, double y, double fwhm);
    void swapX( dComplex *dcData, int lDim, int lAdd );
    bool ParabolicInterpolation(double *xParabola, double *yParabola, double &xMax, double &yMax);
    void fuzzyBinning(double value, double min, double max, int nBins, iPoint2D &iBin, dPoint2D &weightBin );
    double polynomialLegendre(int iPoly, double x);

    void matrixProduct( dMatrix X1, dMatrix X2, dMatrix &X1X2 );
    void matrixProductTranspose1( dMatrix X1, dMatrix X2, dMatrix &X1X2 );
    void matrixProductTranspose2( dMatrix X1, dMatrix X2, dMatrix &X1X2 );
    void matrixPseudoInverse( dMatrix X_12, dMatrix &piX_21 );
    double traceOfSquareMatrix(dMatrix mat);

    void eigen_decomposition(double A[NMAT][NMAT], double V[NMAT][NMAT], double d[NMAT]);
    inline double hypot2(double x, double y) {return qSqrt(x*x+y*y);}
    void tred2(double V[NMAT][NMAT], double d[NMAT], double e[NMAT]);
    void tql2(double V[NMAT][NMAT], double d[NMAT], double e[NMAT]);

    double ibeta(double aa, double bb, double xx);
    double incbcf(double a, double b, double x);
    double incbd(double a, double b, double x);
    double pseries(double a, double b, double x);

    void topDownMergeSort(dVector &array);
    void topDownSplitMerge(int iBegin, int iEnd, dVector &array, dVector &arrayWorking);
    void topDownMerge(int iBegin, int iMiddle, int iEnd, dVector &array, dVector &arrayWorking);
    inline double medianValue(dVector array) {topDownMergeSort(array); return array[array.size()/2];}
    void copyArray(int iBegin, int iEnd, dVector &array, dVector &arrayWorking);

    void topDownMergeSort(dVector &array, iVector &indexArray);
    void topDownSplitMerge(int iBegin, int iEnd, dVector &array, dVector &arrayWorking, iVector &indexArray, iVector &indexWorking);
    void topDownMerge(int iBegin, int iMiddle, int iEnd, dVector &array, dVector &arrayWorking, iVector &indexArray, iVector &indexWorking);
    void copyArray(int iBegin, int iEnd, dVector &array, dVector &arrayWorking, iVector &indexArray, iVector &indexWorking);

    VoxelSet createVoxelList(int iThread, int nThreads, iPoint3D dim);
}
namespace utilString
{
    void errorMessage(QString errorText);
    int decodeSelectionList(QString inputList, bool usePlusOneSyntax, iVector &includeVolume);
    int decodeSelectionList(QString inputList, bool usePlusOneSyntax, bVector &includeVolume);
    QString recodeSelectionList(iVector includeVolume, bool usePlusOneSyntax);
    QString recodeSelectionList(bVector includeVolume, bool usePlusOneSyntax);
    void replaceEnvironmentVariables(QStringList &inputStringList);
    QString replaceEnvironmentVariables(QString inputString);
    QString insertEnvVariables(QString inputString, QStringList templateDirectories);
    QString getFileNameExtension(QString fileName);
    QString getFileNameWithoutExtension(QString fileName);
    QString getFileNameWithoutDirectory(QString fileName);
    QString getFileNameBase(QString fileName);
    QString getDirectoryName(QString fileName);
    bool fileHasExtension(QString fileName);
}
#endif // IO_H
