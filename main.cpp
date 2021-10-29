#include "petmrMain.h"
#include <QApplication>

int main(int argc, char *argv[])
{
    QApplication app(argc, argv);
    MainWindow win;
    win.show();

    const QIcon *cabinetIcon = new QIcon(":/My-Icons/fileCabinet.png");
    app.setWindowIcon(*cabinetIcon);

    return app.exec();
}
