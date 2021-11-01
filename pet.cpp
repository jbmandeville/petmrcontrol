#include <QtWidgets>
#include "petmrMain.h"
#include "ImageIO.h"

void MainWindow::createPETPage()
{
    _petPage = new QWidget();

    /*
    _petRunItemsBox = new QListWidget();
    _petRunItemsBox->setSizePolicy(QSizePolicy::Preferred,QSizePolicy::MinimumExpanding);
    _petRunItemsBox->setSizeAdjustPolicy(QAbstractScrollArea::AdjustToContents);
    connect(_petRunItemsBox, SIGNAL(itemChanged(QListWidgetItem*)),
            this, SLOT(changedfMRIRunCheckBox(QListWidgetItem*)));
    connect(_petRunItemsBox, SIGNAL(itemSelectionChanged()),
            this, SLOT(changedfMRIRunSelection()));
    */

    auto *pageLayout = new QVBoxLayout();
    _petPage->setLayout(pageLayout);
}
