#include <QtMath>
#include <QTime>
#include <QCoreApplication>
#include <QString>
#include <QDebug>
#include <QMessageBox>
#include <QFileDialog>
#include <QFile>
#include "io.h"
#include <fftw3.h>

namespace utilIO
{
    void swap(short *ptr)          { swap_16(ptr);};
    void swap(unsigned short *ptr) { swap_16(ptr);};
    void swap(float *ptr)          { swap_32(ptr);};
    void swap(int  *ptr)           { swap_32(ptr);};

    void swap_16( void *ptr )
    {
        unsigned char uc, *ucptr;

        ucptr = (unsigned char *)ptr;
        uc = *ucptr; *ucptr = *(ucptr+1); *(ucptr+1) = uc;

        return;
    }

    void swap_32( void *ptr )
    {
        unsigned char uc, *ucptr;

        ucptr= (unsigned char*)ptr;
        uc= *ucptr;     *ucptr= *(ucptr+3); *(ucptr+3)= uc;
        uc= *(ucptr+1); *(ucptr+1) = *(ucptr+2); *(ucptr+2)= uc;

        return;
    }

    int machine_byte_order()  // find byte order for cpu
    {
      union
      {
        unsigned char bb[2];
        short         ishort;
      } test;

      int byte_order;

      test.bb[0] = 1;
      test.bb[1] = 0;

      if ( test.ishort == 1 )
        byte_order = 1;
      else
        byte_order = 0;

      return (byte_order);
    }

    int readOverlayFile(VoxelSet &voxelList, QString fileName, overlayColor color, iPoint3D dimSpace)
    {
        QFile file(fileName);
        if (!file.open(QIODevice::ReadOnly | QIODevice::Text))
        {
            qInfo() << "Error reading overlay file " + fileName;
            return(1);
        }
        QTextStream in_stream(&file);

        voxelList.resize(0);
        fVoxel voxel;
        int iLine=1;
        QRegExp rx("[,\\s]");// match a comma or a space
        while (!in_stream.atEnd())
        {
            QString line = in_stream.readLine();
            QString unCommented = line.left(line.indexOf("#"));
            if ( unCommented.isEmpty() )
                continue;
            QStringList valueList = unCommented.split(rx, QString::SkipEmptyParts);
            if ( !valueList.at(0).compare("dimensions",Qt::CaseInsensitive) )
            { // typically this should be on the 1st line, if used
                iPoint3D dim;
                dim.x = valueList.at(1).toInt();
                dim.y = valueList.at(2).toInt();
                dim.z = valueList.at(3).toInt();
                if ( dim.x != dimSpace.x || dim.y != dimSpace.y || dim.z != dimSpace.z )
                {
                    qInfo() << "Error: dimensions of overlay do not match dimensions of image volume.";
                    return(1);
                }
            }
            else
            {
                voxel.x = valueList.at(0).toInt();
                voxel.y = valueList.at(1).toInt();
                voxel.z = valueList.at(2).toInt();
                if ( valueList.size() < 3 )
                    qInfo() << "Error: only read " << valueList.size() << " values on line " << iLine;
                else if ( valueList.size() == 3 )
                    voxel.value = 1.;
                else
                    voxel.value = valueList.at(3).toFloat();
                voxel.color = color;
                voxelList.append(voxel);
            }
            iLine++;
        }
        file.close();
        return(0);
    }

    void delayMS( int millisecondsToWait )
    {
        QTime dieTime = QTime::currentTime().addMSecs( millisecondsToWait );
        while( QTime::currentTime() < dieTime )
        {
            QCoreApplication::processEvents( QEventLoop::AllEvents, 100 );
        }
    }

    QString readTimeTableFile(QString fileName, QStringList &columnNames, dMatrix &table)
    {   // table[time][column] should be allocated previously
        // The time length will not change; the column length may change
        QFile file(fileName);
        if (!file.open(QIODevice::ReadOnly | QIODevice::Text))
        {
            QString errorString = "Error attempting to open file " + fileName;
            return errorString;
        }
        QTextStream in_stream(&file);
        QString line = in_stream.readLine();
        QString unCommented = line.left(line.indexOf("#"));
        QRegExp rx("[,\\s]");// match a comma or a space
        columnNames = unCommented.split(rx, QString::SkipEmptyParts);
        int nColumns = columnNames.size();
        for (int jt=0; jt<table.size(); jt++)
            table[jt].clear();

        int iTime = 0;
        while ( !in_stream.atEnd() )
        {
            QString line = in_stream.readLine();
            QString unCommented = line.left(line.indexOf("#"));
            if ( !unCommented.isEmpty() )
            {
                QStringList valueList = unCommented.split(QRegExp("[,\\s]"), QString::SkipEmptyParts);
                if ( valueList.size() != 0 )
                {
                    if ( valueList.size() < nColumns )
                    {
                        QString errorString = QString("Error: # columns = %1 but expected # = %2 on line %3.").
                                arg(valueList.size()).arg(nColumns).arg(iTime);
                        return errorString;
                    }
                    else
                    {
                        for ( int jColumn=0; jColumn<nColumns; jColumn++)
                        {
                            QString valueString = valueList.at(jColumn);
                            bool ok;
                            double value = valueString.toDouble(&ok);
                            if ( ok )
                                table[iTime].append(value);
                        } // jColumn
                    } // < nColumns
                } // valueSize != 0
            } // !empty
            iTime++;
        } // new line
        file.close();
        return "";
    }

}

////////////////////////////////////////////////////////////////////////////

#ifndef TK_SPLINE_H
#define TK_SPLINE_H

//#include <cstdio>
#include <algorithm>
namespace tk
{

// ---------------------------------------------------------------------
// implementation part, which could be separated into a cpp file
// ---------------------------------------------------------------------


// band_matrix implementation
// -------------------------

void band_matrix::resize(int dim, int n_u, int n_l)
{
    Q_ASSERT(dim>0);
    Q_ASSERT(n_u>=0);
    Q_ASSERT(n_l>=0);
    m_upper.resize(n_u+1);
    m_lower.resize(n_l+1);
    for (int i=0; i<m_upper.size(); i++)
        m_upper[i].resize(dim);
    for (int i=0; i<m_lower.size(); i++)
        m_lower[i].resize(dim);
}

// defines the new operator (), so that we can access the elements
// by A(i,j), index going from i=0,...,dim()-1
double & band_matrix::operator () (int i, int j)
{
    int k=j-i;       // what band is the entry
    Q_ASSERT( (i>=0) && (i<dim()) && (j>=0) && (j<dim()) );
    Q_ASSERT( (-num_lower()<=k) && (k<=num_upper()) );
    // k=0 -> diogonal, k<0 lower left part, k>0 upper right part
    if (k>=0)
        return m_upper[k][i];
    else
        return m_lower[-k][i];
}
double band_matrix::operator () (int i, int j) const
{
    int k=j-i;       // what band is the entry
    Q_ASSERT( (i>=0) && (i<dim()) && (j>=0) && (j<dim()) );
    Q_ASSERT( (-num_lower()<=k) && (k<=num_upper()) );
    // k=0 -> diogonal, k<0 lower left part, k>0 upper right part
    if(k>=0)   return m_upper[k][i];
    else	    return m_lower[-k][i];
}

// LR-Decomposition of a band matrix
void band_matrix::lu_decompose()
{
    int  i_max,j_max;
    int  j_min;
    double x;

    // preconditioning
    // normalize column i so that a_ii=1
    for(int i=0; i<this->dim(); i++) {
        Q_ASSERT(this->operator()(i,i)!=0.0);
        this->saved_diag(i)=1.0/this->operator()(i,i);
        j_min=qMax(0,i-this->num_lower());
        j_max=qMin(this->dim()-1,i+this->num_upper());
        for(int j=j_min; j<=j_max; j++) {
            this->operator()(i,j) *= this->saved_diag(i);
        }
        this->operator()(i,i)=1.0;          // prevents rounding errors
    }

    // Gauss LR-Decomposition
    for(int k=0; k<this->dim(); k++) {
        i_max=qMin(this->dim()-1,k+this->num_lower());  // num_lower not a mistake!
        for(int i=k+1; i<=i_max; i++) {
            Q_ASSERT(this->operator()(k,k)!=0.0);
            x=-this->operator()(i,k)/this->operator()(k,k);
            this->operator()(i,k)=-x;                         // assembly part of L
            j_max=qMin(this->dim()-1,k+this->num_upper());
            for(int j=k+1; j<=j_max; j++) {
                // assembly part of R
                this->operator()(i,j)=this->operator()(i,j)+x*this->operator()(k,j);
            }
        }
    }
}
// solves Ly=b
dVector band_matrix::l_solve(const dVector& b) const
{
    Q_ASSERT( this->dim()==(int)b.size() );
    dVector x(this->dim());
    int j_start;
    double sum;
    for(int i=0; i<this->dim(); i++) {
        sum=0;
        j_start=qMax(0,i-this->num_lower());
        for(int j=j_start; j<i; j++) sum += this->operator()(i,j)*x[j];
        x[i]=(b[i]*this->saved_diag(i)) - sum;
    }
    return x;
}
// solves Rx=y
dVector band_matrix::r_solve(const dVector& b) const
{
    Q_ASSERT( this->dim()==(int)b.size() );
    dVector x(this->dim());
    int j_stop;
    double sum;
    for(int i=this->dim()-1; i>=0; i--) {
        sum=0;
        j_stop=qMin(this->dim()-1,i+this->num_upper());
        for(int j=i+1; j<=j_stop; j++) sum += this->operator()(i,j)*x[j];
        x[i]=( b[i] - sum ) / this->operator()(i,i);
    }
    return x;
}

dVector band_matrix::lu_solve(const dVector& b,
        bool is_lu_decomposed)
{
    Q_ASSERT( this->dim()==(int)b.size() );
    dVector  x,y;
    if(is_lu_decomposed==false) {
        this->lu_decompose();
    }
    y=this->l_solve(b);
    x=this->r_solve(y);
    return x;
}

// spline implementation
// -----------------------

void spline::set_boundary(spline::bd_type left, double left_value,
                          spline::bd_type right, double right_value,
                          bool force_linear_extrapolation)
{
    Q_ASSERT(m_x.size()==0);          // set_points() must not have happened yet
    m_left=left;
    m_right=right;
    m_left_value=left_value;
    m_right_value=right_value;
    m_force_linear_extrapolation=force_linear_extrapolation;
}


void spline::set_points(const dVector& x,
                        const dVector& y, bool cubic_spline)
{
    Q_ASSERT(x.size()==y.size());
    Q_ASSERT(x.size()>2);
    m_x=x;
    m_y=y;
    int   n=x.size();
    // TODO: maybe sort x and y, rather than returning an error
    for(int i=0; i<n-1; i++) {
        Q_ASSERT(m_x[i]<m_x[i+1]);
    }

    if(cubic_spline==true) { // cubic spline interpolation
        // setting up the matrix and right hand side of the equation system
        // for the parameters b[]
        band_matrix A(n,1,1);
        dVector  rhs(n);
        for(int i=1; i<n-1; i++) {
            A(i,i-1)=1.0/3.0*(x[i]-x[i-1]);
            A(i,i)=2.0/3.0*(x[i+1]-x[i-1]);
            A(i,i+1)=1.0/3.0*(x[i+1]-x[i]);
            rhs[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
        }
        // boundary conditions
        if(m_left == spline::second_deriv) {
            // 2*b[0] = f''
            A(0,0)=2.0;
            A(0,1)=0.0;
            rhs[0]=m_left_value;
        } else if(m_left == spline::first_deriv) {
            // c[0] = f', needs to be re-expressed in terms of b:
            // (2b[0]+b[1])(x[1]-x[0]) = 3 ((y[1]-y[0])/(x[1]-x[0]) - f')
            A(0,0)=2.0*(x[1]-x[0]);
            A(0,1)=1.0*(x[1]-x[0]);
            rhs[0]=3.0*((y[1]-y[0])/(x[1]-x[0])-m_left_value);
        } else {
            Q_ASSERT(false);
        }
        if(m_right == spline::second_deriv) {
            // 2*b[n-1] = f''
            A(n-1,n-1)=2.0;
            A(n-1,n-2)=0.0;
            rhs[n-1]=m_right_value;
        } else if(m_right == spline::first_deriv) {
            // c[n-1] = f', needs to be re-expressed in terms of b:
            // (b[n-2]+2b[n-1])(x[n-1]-x[n-2])
            // = 3 (f' - (y[n-1]-y[n-2])/(x[n-1]-x[n-2]))
            A(n-1,n-1)=2.0*(x[n-1]-x[n-2]);
            A(n-1,n-2)=1.0*(x[n-1]-x[n-2]);
            rhs[n-1]=3.0*(m_right_value-(y[n-1]-y[n-2])/(x[n-1]-x[n-2]));
        } else {
            Q_ASSERT(false);
        }

        // solve the equation system to obtain the parameters b[]
        m_b=A.lu_solve(rhs);

        // calculate parameters a[] and c[] based on b[]
        m_a.resize(n);
        m_c.resize(n);
        for(int i=0; i<n-1; i++) {
            m_a[i]=1.0/3.0*(m_b[i+1]-m_b[i])/(x[i+1]-x[i]);
            m_c[i]=(y[i+1]-y[i])/(x[i+1]-x[i])
                   - 1.0/3.0*(2.0*m_b[i]+m_b[i+1])*(x[i+1]-x[i]);
        }
    } else { // linear interpolation
        m_a.resize(n);
        m_b.resize(n);
        m_c.resize(n);
        for(int i=0; i<n-1; i++) {
            m_a[i]=0.0;
            m_b[i]=0.0;
            m_c[i]=(m_y[i+1]-m_y[i])/(m_x[i+1]-m_x[i]);
        }
    }

    // for left extrapolation coefficients
    m_b0 = (m_force_linear_extrapolation==false) ? m_b[0] : 0.0;
    m_c0 = m_c[0];

    // for the right extrapolation coefficients
    // f_{n-1}(x) = b*(x-x_{n-1})^2 + c*(x-x_{n-1}) + y_{n-1}
    double h=x[n-1]-x[n-2];
    // m_b[n-1] is determined by the boundary condition
    m_a[n-1]=0.0;
    m_c[n-1]=3.0*m_a[n-2]*h*h+2.0*m_b[n-2]*h+m_c[n-2];   // = f'_{n-2}(x_{n-1})
    if(m_force_linear_extrapolation==true)
        m_b[n-1]=0.0;
}

double spline::operator() (double x) const
{
    int n=m_x.size();
    // find the closest point m_x[idx] < x, idx=0 even if x<m_x[0]
    dVector::const_iterator it;
    it=std::lower_bound(m_x.begin(),m_x.end(),x);
    int idx=qMax( int(it-m_x.begin())-1, 0);

    double h=x-m_x[idx];
    double interpol;
    if(x<m_x[0]) {
        // extrapolation to the left
        interpol=(m_b0*h + m_c0)*h + m_y[0];
    } else if(x>m_x[n-1]) {
        // extrapolation to the right
        interpol=(m_b[n-1]*h + m_c[n-1])*h + m_y[n-1];
    } else {
        // interpolation
        interpol=((m_a[idx]*h + m_b[idx])*h + m_c[idx])*h + m_y[idx];
    }
    return interpol;
}

} // namespace tk

#endif /* TK_SPLINE_H */

////////////////////////////////////////////////////////////////////////////
namespace utilMath
{

void FFTW1D_Volume( iPoint3D dim, dcVector &volume, int iDimFFT, bool forward )
// iDimFFT: 0 -> x, 1 -> y, 2 -> z
// forward: k-space to image space
{
    int lInner, lOuter1, lOuter2;
    if ( iDimFFT == 0 )
    {
        lInner  = dim.x;
        lOuter1 = dim.y;
        lOuter2 = dim.z;
    }
    else if ( iDimFFT == 1 )
    {
        lOuter1 = dim.x;
        lInner  = dim.y;
        lOuter2 = dim.z;
    }
    else
    {
        lOuter1 = dim.x;
        lOuter2 = dim.y;
        lInner  = dim.z;
    }
    int direction;
    double div;
    if ( forward )
    {
        direction = FFTW_FORWARD;
        div = static_cast<double>(lInner);
    }
    else
    {
        direction = FFTW_BACKWARD;
        div = 1.;
    }

    dComplex *dcData1D   = utilIO::allocate_vector<dComplex>(lInner);
    fftw_complex *fftwData  = static_cast<fftw_complex *>(fftw_malloc(sizeof(fftw_complex)*static_cast<unsigned long>(lInner))); // doubles
    fftw_plan plan_ft1d =fftw_plan_dft_1d(lInner, fftwData, fftwData, direction, FFTW_ESTIMATE);

    for (int j1=0; j1<lOuter1; j1++)
    {
        for (int j2=0; j2<lOuter2; j2++)
        {
            for (int jIn=0; jIn<lInner; jIn++)
            {
                int index;
                if ( iDimFFT == 0 )
                    index = utilIO::index3d(dim, jIn, j1, j2);
                else if ( iDimFFT == 1 )
                    index = utilIO::index3d(dim, j1, jIn, j2);
                else
                    index = utilIO::index3d(dim, j1, j2, jIn);
                dcData1D[jIn].real = volume[index].real;
                dcData1D[jIn].imag = volume[index].imag;
            } // jIn
            // swap on input
            utilMath::swapX( dcData1D, lInner, 0 );
            // Assign the fftw directly (rather than use pointers) in case the byte size differs for the 2 arrays.
            for (int jIn=0; jIn<lInner; jIn++)
            {
                fftwData[jIn][0] = dcData1D[jIn].real;
                fftwData[jIn][1] = dcData1D[jIn].imag;
            }
            // Transform the data.
            fftw_execute(plan_ft1d);
            // Use the original array.
            for (int jIn=0; jIn<lInner; jIn++)
            {
                dcData1D[jIn].real = fftwData[jIn][0];
                dcData1D[jIn].imag = fftwData[jIn][1];
            }
            // The data comes out of the fft routine in wrap-around order.
            utilMath::swapX( dcData1D, lInner, 1 );
            // Put the transformed data back into the original matrix.
            for (int jIn=0; jIn<lInner; jIn++)
            {
                int index;
                if ( iDimFFT == 0 )
                    index = utilIO::index3d(dim, jIn, j1, j2);
                else if ( iDimFFT == 1 )
                    index = utilIO::index3d(dim, j1, jIn, j2);
                else
                    index = utilIO::index3d(dim, j1, j2, jIn);
                volume[index].real = dcData1D[jIn].real / div;
                volume[index].imag = dcData1D[jIn].imag / div;
            } // jIn
        } // j2
    } // j1
    // Free data and plan allocation
    utilIO::delete_vector(dcData1D);
    fftw_destroy_plan(plan_ft1d);
    fftw_free(fftwData);
}

dComplex complexMultiply(dComplex dc1, dComplex dc2, bool star)
{
    // dc1*dc2  or dc1*(dc2*)
    double sign = 1.;
    if ( star ) sign = -1.;

    dComplex dc;
    dc.real = dc1.real * dc2.real - sign * dc1.imag * dc2.imag;
    dc.imag = dc1.real * dc2.imag + sign * dc1.imag * dc2.real;
    return(dc);
}

double computePhase(dComplex dc)
{
    double phase;

    if ( dc.real != 0. )
    {
        phase = qAtan(dc.imag/dc.real);
        if ( dc.real < 0. ) phase += PI;
        if ( phase > PI ) phase -= 2*PI;
    }
    else
    {
        if ( dc.imag == 0. )
            phase = 0.;
        else if ( dc.imag > 0. )
            phase = PI/2.;
        else
            phase = - PI/2.;
    }
    return(phase);
}

double computeAmplitude(dComplex dc)
{
  double amp;
  amp = sqrt(dc.real*dc.real + dc.imag*dc.imag);
  return(amp);
}

double dotProduct(dVector v1, dVector v2)
{
    double sum=0.;
    for (int j=0; j<v1.size(); j++)
        sum += v1[j] * v2[j];
    return sum;
}
double dotProduct(dPoint3D v1, dPoint3D v2)
{
    return v1.x*v2.x + v1.y*v2.y + v1.z*v2.z;
}

bool dInvertSquareMatrix( dMatrix &dSquareMatrix )
{
  //  Purpose
  //  =======
  //  Compute an LU factorization of a square matrix using partial pivoting with row interchanges.
  //
  //  The factorization has the form
  //     dSquareMatrix = P * L * U
  //  where P is a permutation matrix, L is lower triangular with unit
  //  diagonal elements (lower trapezoidal if m > n), and U is upper
  //  triangular (upper trapezoidal if m < n).
  //
  //  This is based upon LAPACK using:
  //  1) dgetf2.f to create the PLU decomposition,
  //  2)
  //
  //  Arguments
  //  =========
  //
  //  dSquareMatrix       (input/output) dimension (nMat,nMat)
  //          On entry, the matrix to be factored.
  //          On exit, the factors L and U from the factorization
  //          dSquareMatrix = P*L*U; the unit diagonal elements of L are not stored.
  //
  //  =====================================================================

    if ( dSquareMatrix.size() != dSquareMatrix[0].size() )
    {
        qWarning() << "Error: only invert square matrix, not " << dSquareMatrix.size() << " by " << dSquareMatrix[0].size() << " matrix.";
        exit(1);
    }
    int nMat = dSquareMatrix.size();
    for (int jRow=0; jRow<nMat; jRow++)
    {
        bool zeroRow = true;
        for (int jCol=0; jCol<nMat; jCol++)
            zeroRow &= (dSquareMatrix[jRow][jCol] == 0.);
        if ( zeroRow ) return false;
    }

    // Allocate working space.
    iVector iPivot;
    dVector dWork;

    iPivot.resize(nMat);
    dWork.resize(nMat);

    ////////////////////// dgetf2.f from LAPACK ////////////////////////////////
    for (int jRow=0; jRow<nMat; jRow++)
    {
        // Choose the largest element to use as the pivot, & test for singularity.
        iPivot[jRow] = jRow;
        double dMax = qAbs(dSquareMatrix[jRow][jRow]);
        for (int j1=jRow+1; j1<nMat; j1++)
        {
              if ( qAbs(dSquareMatrix[j1][jRow]) > dMax )
                {
                  dMax = qAbs(dSquareMatrix[j1][jRow]);
                  iPivot[jRow] = j1;
                }
        }
        if ( dSquareMatrix[iPivot[jRow]][jRow] != 0. )
        {
            // Interchange rows for all columns.
            // E.g., the pivot from the 1st row reorders all columns,
            // the pivot from the 2nd row reorders all columns from the 2nd row down, ...
            if ( iPivot[jRow] != jRow )
            {
              for (int jCol=0; jCol<nMat; jCol++)
              {
                double dTemp                      = dSquareMatrix[jRow][jCol];
                dSquareMatrix[jRow][jCol]         = dSquareMatrix[iPivot[jRow]][jCol];
                dSquareMatrix[iPivot[jRow]][jCol] = dTemp;
              }
            }
            // Compute elements of jRow-th column.
            for (int j1=jRow+1; j1<nMat; j1++)
                dSquareMatrix[j1][jRow] /= dSquareMatrix[jRow][jRow] ;
        }
        else
        {
            qInfo() << "Error: matrix inversion failed on row " << jRow << "with matrix size" << nMat;
            /*
            for (int jRow=0; jRow<nMat; jRow++)
            {
                for (int jCol=0; jCol<nMat; jCol++)
                    qInfo() << "matrix[" << jRow << "][" << jCol << "] =" << dSquareMatrix[jRow][jCol];
            }
            */
            return false; // failed
        }
        // Update trailing submatrix.
        for (int j2=jRow+1; j2<nMat; j2++)
        {
            for (int j1=jRow+1; j1<nMat; j1++)
              dSquareMatrix[j1][j2] -= dSquareMatrix[j1][jRow]*dSquareMatrix[jRow][j2];
        }
    }
    ////////////////////// end of dgetf2.f ////////////////////////////////

    ////////////////////// dtrti2.f from LAPACK ////////////////////////////////
    // Compute the inverse of the upper triangular matrix.
    for (int jCol=0; jCol<nMat; jCol++)
    {
        dSquareMatrix[jCol][jCol] = 1. / dSquareMatrix[jCol][jCol];
        double dDiag = -dSquareMatrix[jCol][jCol];
        // Compute elements 1:j-1 of jCol-th column.
        for ( int jRow=0; jRow<=jCol-1; jRow++)
        {
            for (int j1=0; j1<=jRow-1; j1++)
              dSquareMatrix[j1][jCol] += dSquareMatrix[j1][jRow] * dSquareMatrix[jRow][jCol];
            dSquareMatrix[jRow][jCol] *= dSquareMatrix[jRow][jRow];
        }
        for (int j1=0; j1<=jCol-1; j1++)
            dSquareMatrix[j1][jCol] *= dDiag;
    }
    ////////////////////// end of dtrti2.f ////////////////////////////////

    ////////////////////// dgetri.f from LAPACK ////////////////////////////////
    // Solve the equation inv(dSquareMatrix)*L = inv(U) for inv(dSquareMatrix).
    for (int jCol = nMat-1; jCol>=0; jCol--)
    {
        // Copy current column of L to dWork and replace with zeros.
        for (int jRow=jCol+1; jRow<nMat; jRow++)
        {
          dWork[jRow] = dSquareMatrix[jRow][jCol];
          dSquareMatrix[jRow][jCol] = 0.;
        }
        // Compute current column of inv(dSquareMatrix).
        if (jCol < nMat-1 )
        {
            for ( int jRow=0; jRow<nMat-jCol-1; jRow++ )
            {
                for (int j2=0; j2<nMat; j2++)
                    dSquareMatrix[j2][jCol] -= dSquareMatrix[j2][jCol+jRow+1] * dWork[jCol+jRow+1];;
            }
        }
    }
    // Interchange columns for each row.
    for (int jCol=nMat-2; jCol>=0; jCol--)
    {
        if ( iPivot[jCol] != jCol )
        {
            for (int jRow=0; jRow<nMat; jRow++)
            {
                double dTemp                      = dSquareMatrix[jRow][jCol];
                dSquareMatrix[jRow][jCol]         = dSquareMatrix[jRow][iPivot[jCol]];
                dSquareMatrix[jRow][iPivot[jCol]] = dTemp;
            }
        }
    }
    return true;
    ////////////////////// end of dgetri.f ////////////////////////////////
}

double deltaKernel(double x)
{
    int iX = qRound(x);
    if ( iX == 0 )
        return(1.);
    else
        return(0.);
}

double LanczosKernel(double x, double width)
{ // iX is width of filter (+- iX); good tradeoff for ringing vs. smoothing
    // This filter is almost identical to ham_sinc above
    double lanczos;
    if ( qAbs(x) < 1.0e-5 )    // avoid division by zero
        lanczos = 1.0;
    else if ( qAbs(x/width) > 1. ) // outside range of filter
        lanczos = 0.;
    else
        lanczos = width/SQR(PI*x) * sin(PI*x) * sin(PI*x/width);
    return lanczos;
}

double GaussKernel(double x, double y, double fwhm)
{
  double sigma = 0.42466 * fwhm;
  double arg = - 0.5 * (x*x+y*y) / sigma / sigma;
  double gauss = qExp( arg );
  return gauss;
}

double inverse_Gauss(double x, double y, double fwhm)
{
  double sigma = 0.42466 * fwhm;
  //  arg = - 0.5 * sigma * sigma * (x*x+y*y); //
  double arg = - 2 * PI * PI * sigma * sigma * (x*x+y*y);
  double gauss = qExp( arg );
  return gauss;
}
void swapX( dComplex *dcData, int lDim, int lAdd )
// dcData: double complex data
// lDim:   dimensionality of dcData
// lAdd:   for lDim odd, call as
//         lAdd = 0 before FT
//         lAdd = 1 after FT  (this put array back into original order)
{
    if ( lDim%2 )
    { // the dimensionality is odd
        dComplex *dcStore = utilIO::allocate_vector<dComplex>(lDim);
        for (int j=0; j<lDim; j++)
            dcStore[j] = dcData[j];
        for (int j=0; j<lDim; j++)
        {
            int lShift = j + lDim/2 + lAdd;
            if ( lShift >= lDim )
                lShift -= lDim;
            dcData[lShift] = dcStore[j];
        }
        utilIO::delete_vector(dcStore);
    }
    else
    {  // the dimensionality is even (typical case)
        // The shift becomes a swap.
        for (int j=0; j<lDim/2; j++)
        {
            int lSwap1 = j;
            int lSwap2 = lSwap1 + lDim/2;
            dComplex dc1 = dcData[lSwap1];
            dcData[lSwap1] = dcData[lSwap2];
            dcData[lSwap2] = dc1;
        }
    }
}

double polynomialLegendre(int iPoly, double x)
{ // return Legendre polynomial of x with order iPoly<=9; use x in the range (-1,1)
    double value;
    if ( iPoly == 0 )
        value = 1.;
    else if ( iPoly == 1 )
        value = x;
    else if ( iPoly == 2 )
        value = 0.5 * (3*x*x -1.);
    else if ( iPoly == 3 )
        value = 0.5 * (5.*qPow(x,3) - 3.*x);
    else if ( iPoly == 4 )
        value = 0.125 * (35.*qPow(x,4) - 30.*x*x + 3.);
    else if ( iPoly == 5 )
        value = 0.125 * (63.*qPow(x,5) - 70.*qPow(x,3) + 15.*x);
    else if ( iPoly == 6 )
        value = 1./16. * (231.*qPow(x,6) - 315.*qPow(x,4)+ 105.*x*x - 5.);
    else if ( iPoly == 7 )
        value = 1./16. * (429.*qPow(x,7) - 693.*qPow(x,5) + 315.*qPow(x,3) - 35.);
    else if ( iPoly == 8 )
        value = 1./128. * (6435.*qPow(x,8) - 12012.*qPow(x,6) + 6930.*qPow(x,4) - 1260.*x*x + 35);
    else // if ( iPoly == 9 )
        value = 1./128. * (12155.*qPow(x,9) - 25740.*qPow(x,7) + 18018.*qPow(x,5) - 4620.*qPow(x,3) + 315*x);
    return value;
}

bool ParabolicInterpolation(double *xParabola, double *yParabola, double &xMax, double &yMax)
// return = error state (true=success)
{
    double x0 = xParabola[1];
    double y0 = yParabola[1];
    // Find the maximum of x & y, and then scale each range to +- 1
    double maxX = qMax(qAbs(xParabola[0]-x0),qAbs(xParabola[2]-x0));
    double maxY = qMax(qAbs(yParabola[0]-y0),qAbs(yParabola[2]-y0));
    if ( maxX == 0. )
    {
        FUNC_INFO << xParabola[0] << xParabola[1] << xParabola[2];
        xMax = yMax = 0.;
        return false;
        qFatal("Programming error: maxX = 0 in ParabolicInterpolation()\n");
    }
    else if ( maxY == 0. )
    { // not parabolic: linear instead
        xMax = yMax = 0.;
        return false;
    }
    double xNew[3], yNew[3];
    for (int j1=0; j1<3; j1++)
        xNew[j1] = (xParabola[j1]-x0)/maxX;
    for (int j1=0; j1<3; j1++)
        yNew[j1] = (yParabola[j1]-y0)/maxY;

    // y = a*x^2 + b*x + c
    // ... can be solved for 3 points in matrix form: Y = X * beta
    double xMatrix[3][3];
    xMatrix[0][0] = xNew[0] * xNew[0];
    xMatrix[0][1] = xNew[0];
    xMatrix[0][2] = 1.;
    xMatrix[1][0] = xNew[1] * xNew[1];
    xMatrix[1][1] = xNew[1];
    xMatrix[1][2] = 1.;
    xMatrix[2][0] = xNew[2] * xNew[2];
    xMatrix[2][1] = xNew[2];
    xMatrix[2][2] = 1.;

    // Form the transpose of this matrix
    double xMatrixT[3][3];
    for (int j1=0; j1<3; j1++)
    {
        for (int j2=0; j2<3; j2++)
            xMatrixT[j1][j2] = xMatrix[j2][j1];
    }
    // Calculate XTX
    dMatrix xTx;  xTx.resize(3);
    for (int j1=0; j1<3; j1++)
    {
        xTx[j1].fill(0.,3);
        for (int j2=0; j2<3; j2++)
        {
            xTx[j1][j2] = 0;
            for (int jSum=0; jSum<3; jSum++)
                xTx[j1][j2] += xMatrixT[j1][jSum] * xMatrix[jSum][j2];
        }
    }
    // Invert XTX
    if ( ! utilMath::dInvertSquareMatrix(xTx) )
        return false;
    double xTy[3];
    // Calculate xTy
    for (int j1=0; j1<3; j1++)
    {
        xTy[j1] = 0;
        for (int jsum=0; jsum<3; jsum++)
            xTy[j1] += xMatrixT[j1][jsum] * yNew[jsum];
    }
    // Finally, calculate beta as (XTX)^1 XT*Y
    double beta[3];
    for (int j1=0; j1<3; j1++)
    {
        beta[j1] = 0;
        for (int jsum=0; jsum<3; jsum++)
            beta[j1] += xTx[j1][jsum] * xTy[jsum];
    }
    // y = a*x^2 + b*x + c has a maximum or minimum at x = -b/(2a), unless a = 0
    if ( beta[0] == 0. )
    {
        xMax = yMax = 0.;
        return false;
    }
    else
    {
        xMax = -beta[1]/2./beta[0];
        yMax = beta[0] * xMax * xMax + beta[1] * xMax + beta[2];
        xMax *= maxX;
        yMax *= maxY;
        xMax += x0;
        yMax += y0;
        return true;
    }
}

void fuzzyBinning(double value, double min, double max, int nBins, iPoint2D &iBin, dPoint2D &weightBin )
{ // triangular shaped bins that overlap
    double range  = max - min;
    double offset = value - min;
    double width = range / static_cast<double>(nBins);
    double binValue = offset/width;
    dPoint2D offCenter;
    iBin.x  = qFloor(binValue);  // 0, 1, 2 ... n-1
    offCenter.x = binValue - (iBin.x + 0.5);  // if 0 then center of bin
    weightBin.x = 1. - qAbs(offCenter.x);
    if ( offCenter.x < 0. )
        iBin.y = iBin.x - 1;
    else
        iBin.y = iBin.x + 1;
    offCenter.y = binValue - (iBin.y + 0.5);  // if 0 then center of bin
    weightBin.y = 1. - qAbs(offCenter.y);
    if ( iBin.x < 0 || iBin.x >= nBins ) iBin.x = -1;  // this should never happen
    if ( iBin.y < 0 || iBin.y >= nBins ) iBin.y = -1;
}

void matrixProduct( dMatrix X1, dMatrix X2, dMatrix &X1X2 )
// X1X2_ij = X1_ik * X2_kj
{
    int length1   = X1.size();
    int length2   = X2[0].size();
    int lengthSum = X1[0].size();
    if ( X2.size() != lengthSum )
        qFatal("Inconsistent matrix dimensions in utilMath::matrixProduct");
    X1X2.resize(length1);
    for (int j1=0; j1<X1X2.size(); j1++)
    {
        X1X2[j1].fill(0.,length2);
        for (int j2=0; j2<length2; j2++)
        {
            for (int jSum=0; jSum<lengthSum; jSum++)
                X1X2[j1][j2] += X1[j1][jSum] * X2[jSum][j2];
        }
    }
}

void matrixProductTranspose1( dMatrix X1, dMatrix X2, dMatrix &X1X2 )
// X1X2_ij = (X1_ki)+ * X2_kj
{
    int length1   = X1[0].size();
    int length2   = X2[0].size();
    int lengthSum = X1.size();
    if ( X2.size() != lengthSum )
        qFatal("Inconsistent matrix dimensions in utilMath::matrixProductTranspose1");
    X1X2.resize(length1);
    for (int j1=0; j1<X1X2.size(); j1++)
    {
        X1X2[j1].fill(0.,length2);
        for (int j2=0; j2<length2; j2++)
        {
            for (int jSum=0; jSum<lengthSum; jSum++)
                X1X2[j1][j2] += X1[jSum][j1] * X2[jSum][j2];
        }
    }
}

void matrixProductTranspose2( dMatrix X1, dMatrix X2, dMatrix &X1X2 )
// X1X2_ij = X1_ik * (X2_jk)+
{
    int length1   = X1.size();
    int length2   = X2.size();
    int lengthSum = X1[0].size();
    if ( X2[0].size() != lengthSum )
        qFatal("Inconsistent matrix dimensions in utilMath::matrixProductTranspose2");
    X1X2.resize(length1);
    for (int j1=0; j1<X1X2.size(); j1++)
    {
        X1X2[j1].fill(0.,length2);
        for (int j2=0; j2<length2; j2++)
        {
            for (int jSum=0; jSum<lengthSum; jSum++)
                X1X2[j1][j2] += X1[j1][jSum] * X2[j2][jSum];
        }
    }
}

void matrixPseudoInverse( dMatrix X_12, dMatrix &piX_21 )
{ // Calculate the pseudo-inverse of X_12; note thtat piX_21 does not need to be initialized

  int length2 = X_12[0].size();
  dMatrix XtXm1_22;  XtXm1_22.resize(length2);
  for (int j=0; j<length2; j++)
      XtXm1_22[j].fill(0.,length2);

  // Calculate XTX
  matrixProductTranspose1(X_12, X_12, XtXm1_22);

  // Invert the matrix to get (XTX)^-1
  dInvertSquareMatrix(XtXm1_22);

  // Calculate the pseudo-inverse as (XTX)^-1 XT
  matrixProductTranspose2(XtXm1_22, X_12, piX_21);
}

double traceOfSquareMatrix(dMatrix mat)
{
    if ( mat.size() != mat[0].size() )
        qFatal("Inconsistent matrix dimensions in utilMath::traceOfSquareMatrix");
    double trace=0.;
    for (int j=0; j<mat.size(); j++)
        trace += mat[j][j];
    return(trace);
}

////////////////////////////////////////////////////////////////////////////
//C++ module 'eig3' by Connelly Barnes
//------------------------------------
//
//License: public domain.
//The source files in this directory have been copied from the public domain
//Java matrix library JAMA.  The derived source code is in the public domain as well.
//Eigen decomposition code for symmetric 3x3 matrices, copied from the public domain Java Matrix library JAMA

// Symmetric matrix A => eigenvectors in columns of V, corresponding eigenvalues in d.
// Return evalue[j] in ascending order by magnitude
// Return corresponding evector[i][j] with i=(x,y,z) and j = evalue index
void eigen_decomposition(double A[NMAT][NMAT], double V[NMAT][NMAT], double d[NMAT])
{
    double e[NMAT];
    int i, j;
    for (i = 0; i < NMAT; i++)
    {
        for (j = 0; j < NMAT; j++)
        {
            V[i][j] = A[i][j];
        }
    }
    tred2(V, d, e);
    tql2(V, d, e);
}

// Symmetric Householder reduction to tridiagonal form.
void tred2(double V[NMAT][NMAT], double d[NMAT], double e[NMAT])
{
    //  This is derived from the Algol procedures tred2 by
    //  Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
    //  Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
    //  Fortran subroutine in EISPACK.
    int i, j, k;

    for (j = 0; j < NMAT; j++)
    {
        d[j] = V[NMAT-1][j];
    }

    // Householder reduction to tridiagonal form.

    for (i = NMAT-1; i > 0; i--)
    {
        // Scale to avoid under/overflow.

        double scale = 0.0;
        double h = 0.0;
        for (k = 0; k < i; k++)
        {
            scale = scale + qAbs(d[k]);
        }
        if (scale == 0.0)
        {
            e[i] = d[i-1];
            for (j = 0; j < i; j++)
            {
                d[j] = V[i-1][j];
                V[i][j] = 0.0;
                V[j][i] = 0.0;
            }
        }
        else
        {

            // Generate Householder vector.
            for (k = 0; k < i; k++)
            {
                d[k] /= scale;
                h += d[k] * d[k];
            }
            double f = d[i-1];
            double g = qSqrt(h);
            if (f > 0) {
                g = -g;
            }
            e[i] = scale * g;
            h = h - f * g;
            d[i-1] = f - g;
            for (j = 0; j < i; j++)
            {
                e[j] = 0.0;
            }

            // Apply similarity transformation to remaining columns.

            for (j = 0; j < i; j++)
            {
                f = d[j];
                V[j][i] = f;
                g = e[j] + V[j][j] * f;
                for (k = j+1; k <= i-1; k++)
                {
                    g += V[k][j] * d[k];
                    e[k] += V[k][j] * f;
                }
                e[j] = g;
            }
            f = 0.0;
            for (j = 0; j < i; j++)
            {
                e[j] /= h;
                f += e[j] * d[j];
            }
            double hh = f / (h + h);
            for (j = 0; j < i; j++)
            {
                e[j] -= hh * d[j];
            }
            for (j = 0; j < i; j++)
            {
                f = d[j];
                g = e[j];
                for (k = j; k <= i-1; k++)
                {
                    V[k][j] -= (f * e[k] + g * d[k]);
                }
                d[j] = V[i-1][j];
                V[i][j] = 0.0;
            }
        }
        d[i] = h;
    }

    // Accumulate transformations.

    for (i = 0; i < NMAT-1; i++)
    {
        V[NMAT-1][i] = V[i][i];
        V[i][i] = 1.0;
        double h = d[i+1];
        if (h != 0.0) {
            for (k = 0; k <= i; k++)
            {
                d[k] = V[k][i+1] / h;
            }
            for (j = 0; j <= i; j++)
            {
                double g = 0.0;
                for (k = 0; k <= i; k++)
                {
                    g += V[k][i+1] * V[k][j];
                }
                for (k = 0; k <= i; k++) {
                    V[k][j] -= g * d[k];
                }
            }
        }
        for (k = 0; k <= i; k++)
        {
            V[k][i+1] = 0.0;
        }
    }
    for (j = 0; j < NMAT; j++)
    {
        d[j] = V[NMAT-1][j];
        V[NMAT-1][j] = 0.0;
    }
    V[NMAT-1][NMAT-1] = 1.0;
    e[0] = 0.0;
}

// Symmetric tridiagonal QL algorithm.
void tql2(double V[NMAT][NMAT], double d[NMAT], double e[NMAT])
{
    //  This is derived from the Algol procedures tql2, by
    //  Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
    //  Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
    //  Fortran subroutine in EISPACK.

    int i, j, k, l;

    for (i = 1; i < NMAT; i++)
    {
        e[i-1] = e[i];
    }
    e[NMAT-1] = 0.0;

    double f = 0.0;
    double tst1 = 0.0;
    double eps = qPow(2.0,-52.0);
    for (l = 0; l < NMAT; l++) {

        // Find small subdiagonal element

        tst1 = qMax(tst1,qAbs(d[l]) + qAbs(e[l]));
        int m = l;
        while (m < NMAT)
        {
            if (qAbs(e[m]) <= eps*tst1) {
                break;
            }
            m++;
        }

        // If m == l, d[l] is an eigenvalue,
        // otherwise, iterate.

        if (m > l)
        {
            int iter = 0;
            do {
                iter = iter + 1;  // (Could check iteration count here.)

                // Compute implicit shift

                double g = d[l];
                double p = (d[l+1] - g) / (2.0 * e[l]);
                double r = hypot2(p,1.0);
                if (p < 0)
                {
                    r = -r;
                }
                d[l] = e[l] / (p + r);
                d[l+1] = e[l] * (p + r);
                double dl1 = d[l+1];
                double h = g - d[l];
                for (i = l+2; i < NMAT; i++)
                {
                    d[i] -= h;
                }
                f = f + h;

                // Implicit QL transformation.

                p = d[m];
                double c = 1.0;
                double c2 = c;
                double c3 = c;
                double el1 = e[l+1];
                double s = 0.0;
                double s2 = 0.0;
                for (i = m-1; i >= l; i--)
                {
                    c3 = c2;
                    c2 = c;
                    s2 = s;
                    g = c * e[i];
                    h = c * p;
                    r = hypot2(p,e[i]);
                    e[i+1] = s * r;
                    s = e[i] / r;
                    c = p / r;
                    p = c * d[i] - s * g;
                    d[i+1] = h + s * (c * g + s * d[i]);

                    // Accumulate transformation.

                    for (k = 0; k < NMAT; k++)
                    {
                        h = V[k][i+1];
                        V[k][i+1] = s * V[k][i] + c * h;
                        V[k][i] = c * V[k][i] - s * h;
                    }
                }
                p = -s * s2 * c3 * el1 * e[l] / dl1;
                e[l] = s * p;
                d[l] = c * p;

                // Check for convergence.

            } while (qAbs(e[l]) > eps*tst1);
        }
        d[l] = d[l] + f;
        e[l] = 0.0;
    }

    // Sort eigenvalues and corresponding vectors.

    for (i = 0; i < NMAT-1; i++)
    {
        int k = i;
        double p = d[i];
        for (j = i+1; j < NMAT; j++)
        {
            if (d[j] < p) {
                k = j;
                p = d[j];
            }
        }
        if (k != i)
        {
            d[k] = d[i];
            d[i] = p;
            for (j = 0; j < NMAT; j++)
            {
                p = V[j][i];
                V[j][i] = V[j][k];
                V[j][k] = p;
            }
        }
    }
}


// The following 4 c-code functions:
//                          1) ibeta
//                          2) incbcf
//                          3) incbd
//                          4) pseries
// are modified versions (derivative code) of public domain software.
// This software distribution was originally based upon the
// Cephes Mathematical Library from www.netlib.org/cephes. The particular
// code underlying this derivative version was obtained from
//         www.codeproject.com/KB/cs/SpecialFunction.aspx?display=PRINT
// which uses the "The Code Project Open License" at
//                www.codeproject.com/info/cpol10.aspx
// The license states that ...
// " The main points subject to the terms of the License are:
//  * Source Code and Executable Files can be used in commercial applications;
//  * Source Code and Executable Files can be redistributed; and
//  * Source Code can be modified to create derivative works.
//  * No claim of suitability, guarantee, or any warranty whatsoever is provided.
//  The software is provided "as-is".
//  * The Article(s) accompanying the Work may not be distributed or republished without the Author's consent
//
// The header to the original c code is reproduced below.
////////////////////////////////////////////////////////////////////////////
/*
**************************************************************************
**
**    Class  SpecialFunction (C#)
**
**************************************************************************
**    Copyright (C) 1984 Stephen L. Moshier (original C version - Cephes Math Library)
**    Copyright (C) 1996 Leigh Brookshaw	(Java version)
**    Copyright (C) 2005 Miroslav Stampar	(C# version [->this<-])
**
**    This program is free software; you can redistribute it and/or modify
**    it under the terms of the GNU General Public License as published by
**    the Free Software Foundation; either version 2 of the License, or
**    (at your option) any later version.
**
**    This program is distributed in the hope that it will be useful,
**    but WITHOUT ANY WARRANTY; without even the implied warranty of
**    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
**    GNU General Public License for more details.
**
**    You should have received a copy of the GNU General Public License
**    along with this program; if not, write to the Free Software
**    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
**************************************************************************
**
**    This class is an extension of System.Math. It includes a number
**    of special functions not found in the Math class.
**
*************************************************************************
** This class contains physical constants and special functions not found
** in the System.Math class.
** Like the System.Math class this class is final and cannot be
** subclassed.
** All physical constants are in cgs units.
** NOTE: These special functions do not necessarily use the fastest
** or most accurate algorithms.
**
** @version $Revision: 1.8 $, $Date: 2005/09/12 09:52:34 $
**/

#define MACHEP 1.11022302462515654042E-16
#define MAXLOG 7.09782712893383996732E2
#define MINLOG -7.451332191019412076235E2
#define MAXGAM 171.624376956302725

//prototypes
/*
double incbcf(double a, double b, double x);
double incbd(double a, double b, double x);
double pseries(double a, double b, double x);
double tgamma(double x);
*/

double ibeta(double aa, double bb, double xx)
{
    double a, b, t, x, xc, w, y;
    bool flag;

    if (aa <= 0.0 || bb <= 0.0)
    {
        printf("ibeta: Domain error!, a = %g, b = %g\n",aa,bb);
        exit(1);
    }
    if ((xx <= 0.0) || (xx >= 1.0))
    {
        if (xx == 0.0) return 0.0;
        if (xx == 1.0) return 1.0;
    }
    flag = false;
    if ((bb*xx) <= 1.0 && xx <= 0.95)
    {
        t = pseries(aa, bb, xx);
        return t;
    }

    w = 1.0 - xx;

    /* Reverse a and b if x is greater than the mean. */
    if (xx > (aa/(aa + bb)))
    {
        flag = true;
        a = bb;
        b = aa;
        xc = xx;
        x = w;
    }
    else
    {
        a = aa;
        b = bb;
        xc = w;
        x = xx;
    }

    if (flag && (b*x) <= 1.0 && x <= 0.95)
    {
        t = pseries(a, b, x);
        if (t <= MACHEP)
            t = 1.0 - MACHEP;
        else
            t = 1.0 - t;
        return t;
    }

    /* Choose expansion for better convergence. */
    y = x*(a + b - 2.0) - (a - 1.0);
    if (y < 0.0)
        w = incbcf(a, b, x);
    else
        w = incbd(a, b, x)/xc;

    /* Multiply w by the factor
     a      b   _             _     _
     x  (1-x)   | (a+b) / ( a | (a) | (b) ) .   */

    y = a*qLn(x);
    t = b*qLn(xc);
    if ((a + b) < MAXGAM && qAbs(y) < MAXLOG && qAbs(t) < MAXLOG)
    {
        t = qPow(xc, b);
        t *= qPow(x, a);
        t /= a;
        t *= w;
        t *= std::tgamma(a + b)/(std::tgamma(a)*std::tgamma(b));
        if (flag)
        {
            if (t <= MACHEP) t = 1.0 - MACHEP;
            else t = 1.0 - t;
        }
        return t;
    }
    /* Resort to logarithms.  */
    y += t + std::lgamma(a + b) - std::lgamma(a) - std::lgamma(b);
    y += qLn(w/a);
    if (y < MINLOG)
        t = 0.0;
    else
        t = qExp(y);

    if (flag)
    {
        if (t <= MACHEP) t = 1.0 - MACHEP;
        else t = 1.0 - t;
    }
    return t;
}

/// Returns the continued fraction expansion #1 for incomplete beta integral.
double incbcf(double a, double b, double x)
{
    double xk, pk, pkm1, pkm2, qk, qkm1, qkm2;
    double k1, k2, k3, k4, k5, k6, k7, k8;
    double r, t, ans, thresh;
    int n;
    double big = 4.503599627370496e15;
    double biginv = 2.22044604925031308085e-16;

    k1 = a;
    k2 = a + b;
    k3 = a;
    k4 = a + 1.0;
    k5 = 1.0;
    k6 = b - 1.0;
    k7 = k4;
    k8 = a + 2.0;

    pkm2 = 0.0;
    qkm2 = 1.0;
    pkm1 = 1.0;
    qkm1 = 1.0;
    ans = 1.0;
    r = 1.0;
    n = 0;
    thresh = 3.0*MACHEP;
    do
    {
        xk = -(x*k1*k2)/(k3*k4);
        pk = pkm1 + pkm2*xk;
        qk = qkm1 + qkm2*xk;
        pkm2 = pkm1;
        pkm1 = pk;
        qkm2 = qkm1;
        qkm1 = qk;

        xk = (x*k5*k6)/(k7*k8);
        pk = pkm1 + pkm2*xk;
        qk = qkm1 + qkm2*xk;
        pkm2 = pkm1;
        pkm1 = pk;
        qkm2 = qkm1;
        qkm1 = qk;

        if (qk != 0.) r = pk/qk;
        if (r != 0.)
        {
            t = qAbs((ans - r)/r);
            ans = r;
        }
        else
            t = 1.0;

        if (t < thresh) return ans;

        k1 += 1.0;
        k2 += 1.0;
        k3 += 2.0;
        k4 += 2.0;
        k5 += 1.0;
        k6 -= 1.0;
        k7 += 2.0;
        k8 += 2.0;

        if ((qAbs(qk) + qAbs(pk)) > big)
        {
            pkm2 *= biginv;
            pkm1 *= biginv;
            qkm2 *= biginv;
            qkm1 *= biginv;
        }
        if ((qAbs(qk) < biginv) || (qAbs(pk) < biginv))
        {
            pkm2 *= big;
            pkm1 *= big;
            qkm2 *= big;
            qkm1 *= big;
        }
    } while (++n < 300);

    return ans;
}

/// Returns the continued fraction expansion #2 for incomplete beta integral.
double incbd(double a, double b, double x)
{
    double xk, pk, pkm1, pkm2, qk, qkm1, qkm2;
    double k1, k2, k3, k4, k5, k6, k7, k8;
    double r, t, ans, z, thresh;
    int n;
    double big = 4.503599627370496e15;
    double biginv = 2.22044604925031308085e-16;

    k1 = a;
    k2 = b - 1.0;
    k3 = a;
    k4 = a + 1.0;
    k5 = 1.0;
    k6 = a + b;
    k7 = a + 1.0;
    ;
    k8 = a + 2.0;

    pkm2 = 0.0;
    qkm2 = 1.0;
    pkm1 = 1.0;
    qkm1 = 1.0;
    z = x/(1.0 - x);
    ans = 1.0;
    r = 1.0;
    n = 0;
    thresh = 3.0*MACHEP;
    do
    {
        xk = -(z*k1*k2)/(k3*k4);
        pk = pkm1 + pkm2*xk;
        qk = qkm1 + qkm2*xk;
        pkm2 = pkm1;
        pkm1 = pk;
        qkm2 = qkm1;
        qkm1 = qk;

        xk = (z*k5*k6)/(k7*k8);
        pk = pkm1 + pkm2*xk;
        qk = qkm1 + qkm2*xk;
        pkm2 = pkm1;
        pkm1 = pk;
        qkm2 = qkm1;
        qkm1 = qk;

        if (qk != 0.) r = pk/qk;
        if (r != 0.)
        {
            t = qAbs((ans - r)/r);
            ans = r;
        }
        else
            t = 1.0;

        if (t < thresh) return ans;

        k1 += 1.0;
        k2 -= 1.0;
        k3 += 2.0;
        k4 += 2.0;
        k5 += 1.0;
        k6 += 1.0;
        k7 += 2.0;
        k8 += 2.0;

        if ((qAbs(qk) + qAbs(pk)) > big)
        {
            pkm2 *= biginv;
            pkm1 *= biginv;
            qkm2 *= biginv;
            qkm1 *= biginv;
        }
        if ((qAbs(qk) < biginv) || (qAbs(pk) < biginv))
        {
            pkm2 *= big;
            pkm1 *= big;
            qkm2 *= big;
            qkm1 *= big;
        }
    } while (++n < 300);

    return ans;
}

/// Returns the power series for incomplete beta integral. Use when b*x is small and x not too close to 1.
double pseries(double a, double b, double x)
{
    double s, t, u, v, n, t1, z, ai;

    ai = 1.0/a;
    u = (1.0 - b)*x;
    v = u/(a + 1.0);
    t1 = v;
    t = u;
    n = 2.0;
    s = 0.0;
    z = MACHEP*ai;
    while (qAbs(v) > z)
    {
        u = (n - b)*x/n;
        t *= u;
        v = t/(a + n);
        s += v;
        n += 1.0;
    }
    s += t1;
    s += ai;

    u = a*qLn(x);
    if ((a + b) < MAXGAM && qAbs(u) < MAXLOG)
    {
        t = std::tgamma(a + b)/(std::tgamma(a)*std::tgamma(b));
        s = s*t*qPow(x, a);
    }
    else
    {
        t = std::lgamma(a + b) - std::lgamma(a) - std::lgamma(b) + u + qLn(s);
        if (t < MINLOG) s = 0.0;
        else s = qExp(t);
    }
    return s;
}

void topDownMergeSort(dVector &array)
{ // array has the items to sort
    dVector arrayWorking;   arrayWorking.resize(array.size());
    topDownSplitMerge(0, array.size(), array, arrayWorking);
}
void topDownSplitMerge(int iBegin, int iEnd, dVector &array, dVector &arrayWorking)
{
    int iMiddle;
    if (iEnd - iBegin < 2)                      // if run size == 1
        return;                                 //   consider it sorted
    // recursively split runs into two halves until run size == 1,
    // then merge them and return back up the call chain
    iMiddle = (iEnd + iBegin) / 2;              // iMiddle = mid point
    topDownSplitMerge(iBegin,  iMiddle, array, arrayWorking);  // split / merge left  half
    topDownSplitMerge(iMiddle,    iEnd, array, arrayWorking);  // split / merge right half
    topDownMerge(iBegin, iMiddle, iEnd, array, arrayWorking);  // merge the two half runs
    copyArray(iBegin, iEnd, array, arrayWorking);              // copy the merged runs back to array
}
//  left half is array[iBegin :iMiddle-1]
// right half is array[iMiddle:iEnd-1   ]
void topDownMerge(int iBegin, int iMiddle, int iEnd,
                  dVector &array, dVector &arrayWorking)
{
    int i0, i1, j;
    i0 = iBegin; i1 = iMiddle;
    // While there are elements in the left or right runs
    for (j = iBegin; j < iEnd; j++)
    {
        // If left run head exists and is <= existing right run head.
        if (i0 < iMiddle && (i1 >= iEnd || array[i0] <= array[i1]))
        {
            arrayWorking[j] = array[i0];
            i0 = i0 + 1;
        }
        else
        {
            arrayWorking[j] = array[i1];
            i1 = i1 + 1;
        }
    }
}
void copyArray(int iBegin, int iEnd, dVector &array, dVector &arrayWorking)
{
    int j;
    for(j = iBegin; j < iEnd; j++)
        array[j] = arrayWorking[j];
}

void topDownMergeSort(dVector &array, iVector &indexArray)
{   // array has the items to sort; indexArray carries along the sorted indices for application to other vectors
    // initialize indexArray: allocate and fill
    indexArray.resize(array.size());
    for (int j=0; j<array.size(); j++)
        indexArray[j] = j;

    dVector arrayWorking;   arrayWorking.resize(array.size());
    iVector indexWorking;   indexWorking.resize(indexArray.size());
    topDownSplitMerge(0, array.size(), array, arrayWorking, indexArray, indexWorking);
}
// iBegin is inclusive; iEnd is exclusive (array[iEnd] is not in the set)
void topDownSplitMerge(int iBegin, int iEnd, dVector &array, dVector &arrayWorking, iVector &indexArray, iVector &indexWorking)
{
    int iMiddle;
    if (iEnd - iBegin < 2)                      // if run size == 1
        return;                                 //   consider it sorted
    // recursively split runs into two halves until run size == 1,
    // then merge them and return back up the call chain
    iMiddle = (iEnd + iBegin) / 2;              // iMiddle = mid point
    topDownSplitMerge(iBegin,  iMiddle, array, arrayWorking, indexArray, indexWorking);  // split / merge left  half
    topDownSplitMerge(iMiddle,    iEnd, array, arrayWorking, indexArray, indexWorking);  // split / merge right half
    topDownMerge(iBegin, iMiddle, iEnd, array, arrayWorking, indexArray, indexWorking);  // merge the two half runs
    copyArray(iBegin, iEnd, array, arrayWorking, indexArray, indexWorking);              // copy the merged runs back to array
}
//  left half is array[iBegin :iMiddle-1]
// right half is array[iMiddle:iEnd-1   ]
void topDownMerge(int iBegin, int iMiddle, int iEnd,
                             dVector &array, dVector &arrayWorking,
                             iVector &indexArray, iVector &indexWorking)
{
    int i0, i1, j;
    i0 = iBegin; i1 = iMiddle;
    // While there are elements in the left or right runs
    for (j = iBegin; j < iEnd; j++)
    {
        // If left run head exists and is <= existing right run head.
        if (i0 < iMiddle && (i1 >= iEnd || array[i0] <= array[i1]))
        {
            arrayWorking[j] = array[i0];
            indexWorking[j] = indexArray[i0];
            i0 = i0 + 1;
        }
        else
        {
            arrayWorking[j] = array[i1];
            indexWorking[j] = indexArray[i1];
            i1 = i1 + 1;
        }
    }
}
void copyArray(int iBegin, int iEnd, dVector &array, dVector &arrayWorking, iVector &indexArray, iVector &indexWorking)
{
    int j;
    for(j = iBegin; j < iEnd; j++)
    {
        array[j]       = arrayWorking[j];
        indexArray[j] = indexWorking[j];
    }
}

VoxelSet createVoxelList(int iThread, int nThreads, iPoint3D dim)
{
    VoxelSet voxelList;  voxelList.resize(0);
    int iVoxel=0;
    fVoxel voxel;  voxel.value = 0.;  voxel.color = Color_Undefined;
    for ( voxel.z=0; voxel.z<dim.z; voxel.z++)
    {
        for ( voxel.y=0; voxel.y<dim.y; voxel.y++)
        {
            for ( voxel.x=0; voxel.x<dim.x; voxel.x++)
            {
                if ( (iVoxel-iThread)%nThreads == 0 )
                    voxelList.append(voxel);
                iVoxel++;
            }
        }
    }
    return voxelList;
}

} // namespace utilMath
////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////
namespace utilString
{
void errorMessage(QString errorText)
{
    QMessageBox msgBox;  msgBox.setIcon(QMessageBox::Critical);
    msgBox.setText(errorText);
    msgBox.exec();
    qInfo() << errorText;
}

int decodeSelectionList(QString inputList, bool usePlusOneSyntax, iVector &includeVolume)
{ // includeVolume should be previously allocated to the desired dimension
//    FUNC_ENTER << usePlusOneSyntax << includeVolume;
    bVector includeBool;
    for ( int jt=0; jt<includeVolume.size(); jt++ )
    {
        if ( includeVolume[jt] == 0 )
            includeBool.append(false);
        else
            includeBool.append(true);
    }
    int error = decodeSelectionList(inputList, usePlusOneSyntax, includeBool);
    if ( error ) return error;
    for ( int jt=0; jt<includeBool.size(); jt++ )
    {
        if ( includeBool[jt] )
            includeVolume.replace(jt,1);
        else
            includeVolume.replace(jt,0);
    }
    return 0;
}

int decodeSelectionList(QString inputList, bool usePlusOneSyntax, bVector &includeVolume)
{ // includeVolume should be previously defined with the desired dimension; values with be binary for inclusion
    FUNC_ENTER << usePlusOneSyntax << includeVolume;
    QStringList stringList = inputList.split(QRegExp("[ ]"), QString::SkipEmptyParts);
    if ( inputList.isEmpty() || stringList.count() == 0 ) return 0;

    inputList.replace(QRegExp(" "), "-");
    stringList = inputList.split(QRegExp("[,]"), QString::SkipEmptyParts);

    int iOffset = 0;
    if ( usePlusOneSyntax ) iOffset = -1;
    // Each string in the list should have values separated by dashes or spaces
    for (int jString=0; jString<stringList.count(); jString++)
    {
        QString dashString = stringList.at(jString);
        QStringList valueList = dashString.split(QRegExp("[-]"), QString::SkipEmptyParts);
        if ( valueList.count() == 1 )
            valueList.append(valueList.at(0));
        else if ( valueList.count() != 2 )
            return 1;
        bool ok1, ok2;
        int iLow = valueList.at(0).toInt(&ok1) + iOffset;
        int iHigh = valueList.at(1).toInt(&ok2) + iOffset;
        if ( valueList.at(1) == "n" )
        {
            iHigh = includeVolume.size()-1;
            ok2 = true;
        }
        if ( ok1 && ok2 )
        {
            for (int j=iLow; j<=iHigh; j++)
            {
                if ( j >= 0 && j < includeVolume.size() )
                    includeVolume[j] = true;
                else
                    return 1;
            }
        }
        else
            return 1;
    }
    /*
    if ( !includeVolume.contains(true) )  // must contain at least 1 true value
        return 1;
    else
        return 0;
        */
    return 0;
}

QString recodeSelectionList(iVector includeVolume, bool usePlusOneSyntax)
{
    bVector includeBool;
    for ( int jt=0; jt<includeVolume.size(); jt++ )
    {
        if ( includeVolume[jt] == 0 )
            includeBool.append(false);
        else
            includeBool.append(true);
    }
    return recodeSelectionList(includeBool, usePlusOneSyntax);
}

QString recodeSelectionList(bVector includeVolume, bool usePlusOneSyntax)
{
    QString selectionString = "";
    int it=0;
    int iOffset = 0;
    if ( usePlusOneSyntax ) iOffset = +1;
    while ( it<includeVolume.size() )
    {
        if ( includeVolume[it] )
        {
            int iLow  = it + iOffset;
            while ( it<includeVolume.size() && includeVolume[it] )
                it++;
            int iHigh = it-1 + iOffset;
            if ( iLow == iHigh )
            {
                QString oneValue;  oneValue.setNum(iLow);
                if ( selectionString.length() == 0 )
                    selectionString.append(oneValue);
                else
                    selectionString.append(","+oneValue);
            }
            else
            {
                QString lowString;  lowString.setNum(iLow);
                QString highString; highString.setNum(iHigh);
                if ( selectionString.length() == 0 )
                    selectionString.append(lowString+"-"+highString);
                else
                    selectionString.append(","+lowString+"-"+highString);
            }
        }
        it++;
    }
    return selectionString;
}

void replaceEnvironmentVariables(QStringList &inputStringList)
{ // Replace environment variables in all strings in the list
    for ( int jString=0; jString<inputStringList.count(); jString++)
    {
        QString inputString = inputStringList.at(jString);
        inputStringList[jString] = replaceEnvironmentVariables(inputString);
    }
}

QString replaceEnvironmentVariables(QString inputString)
{ // Return a new string with the environment variables replaced. Do not alter the input string.
    QString outputString = inputString;

    int iDollar = outputString.indexOf('$');
    if ( iDollar < 0 )
        return outputString;

    QStringList envList(QProcess::systemEnvironment());
    FUNC_INFO << envList;
    int endoutputString = outputString.length()-1;
    while ( (iDollar = outputString.indexOf('$', iDollar)) != -1 )
    {
        int startShellVariable  = iDollar;
        QString shellVariable = outputString.right(endoutputString-startShellVariable);
        bool allowed = true;
        int iCount=0;
        while ( iCount < shellVariable.length() && allowed )
        {
            int iAscii = shellVariable.at(iCount).toLatin1();
            allowed = iAscii > 64 && iAscii < 91;   // capital letters
            allowed |= iAscii > 96 && iAscii < 123; // lower-case letter
            allowed |= iAscii > 47 && iAscii < 58;  // numbers
            allowed |= iAscii == 95;                // underscore
            if ( allowed ) iCount++;
        }
        shellVariable = shellVariable.left(iCount);
        // Given a potential shell variable (based upon allowed characters), test all env variables
        bool found = false;
        for ( int jString=0; jString<envList.count(); jString++)
        {
            QString envString = envList.at(jString);
            int equal = envString.indexOf('=');
            int end   = envString.length()-1;
            QString envVariable = envString.left(equal);
            QString envValue    = envString.right(end-equal);
            if ( shellVariable.compare(envVariable) == 0 )
            {
                found = true;
                outputString.replace(iDollar,iCount+1,envValue);  // +1 for $ sign
            }
            while ( outputString.contains("//") )
            {
                int iStart = outputString.indexOf("//");
                outputString.replace(iStart,2,'/');
            }
        }
        if ( ! found )
        {
            qInfo() << "Error: variable" << shellVariable << "not found in shell environment.";
            exit(1);
        }
        ++iDollar;
    }
    return outputString;
}

QString insertEnvVariables(QString inputString, QStringList templateDirectories)
{ // Return a new string with the environment variables replaced. Do not alter the input string.
    FUNC_ENTER << inputString;
    QString outputString = inputString;

    int nTemplates = templateDirectories.count() / 2;
    for (int jTemplate=0; jTemplate<nTemplates; jTemplate++)
    {
        QString envVariable = templateDirectories.at(2*jTemplate);
        QString envValue    = templateDirectories.at(2*jTemplate+1);
        FUNC_INFO << "envVariable" << envVariable << "envValue" << envValue;
        if ( outputString.contains(envValue) )
        {
            FUNC_INFO << "contains!!";
            QString envVariableWithDollarAndSlash = QString("$%1/").arg(envVariable);
            outputString.replace(envValue,envVariableWithDollarAndSlash);
        }
    }
    return outputString;
}

bool fileHasExtension(QString fileName)
{
    // A dot could take these forms:
    // "file.dat"
    // "../file.dat"
    // "../file
    //
    // Determine where a slash exists (end of string if missing)
    int dot   = fileName.lastIndexOf(".");
    int slash = fileName.lastIndexOf("/");
    // Determine the file extension
    bool extension = false;
    if ( (dot>=0 && slash==-1) || (dot>=0 && dot>slash) ) // The file has an extension.
        extension = true;
    return extension;
}
QString getFileNameExtension(QString fileName)
{
    // A dot could take these forms:
    // "file.dat"
    // "../file.dat"
    // "../file
    //
    // Determine where a slash exists (end of string if missing)
    int dot   = fileName.lastIndexOf(".");
    int slash = fileName.lastIndexOf("/");
    int end   = fileName.length()-1;
    // Determine the file extension
    QString extension = "";
    if ( (dot>=0 && slash==-1) || (dot>=0 && dot>slash) ) // The file has an extension.
        extension = fileName.right(end-dot);
    return extension;
}
QString getFileNameWithoutExtension(QString fileName)
{
    // A dot could take these forms:
    // "file.dat"
    // "../file.dat"
    // "../file
    //
    // Determine where a slash exists (end of string if missing)
    int dot   = fileName.lastIndexOf(".");
    int slash = fileName.lastIndexOf("/");
    // Determine the file extension
    QString extensionlessName = fileName;
    if ( (dot>=0 && slash==-1) || (dot>=0 && dot>slash) ) // The file has an extension.
        extensionlessName = fileName.left(dot);
    return extensionlessName;
}
QString getFileNameWithoutDirectory(QString fileName)
{
    QString shortFileName = fileName;
    int slash = shortFileName.lastIndexOf("/");
    if ( slash >= 0 )
    {
        int end   = shortFileName.length()-1;
        shortFileName = shortFileName.right(end-slash);
    }
    return shortFileName;
}
QString getFileNameBase(QString fileName)
{
    QString baseName = getFileNameWithoutExtension(fileName);
    int slash = baseName.lastIndexOf("/");
    if ( slash >= 0 )
    {
        int end   = baseName.length()-1;
        baseName = baseName.right(end-slash);
    }
    return baseName;
}
QString getDirectoryName(QString fileName)
{
    QString directoryName = fileName;
    int slash = directoryName.lastIndexOf("/");
    if ( slash >= 0 )
        directoryName = directoryName.left(slash);
    else
        directoryName = ".";
    return directoryName;
}

} // namespace utilString
////////////////////////////////////////////////////////////////////////////
