#include "mainwindow.h"

mainwindow::mainwindow(QWidget *parent)
    : QMainWindow(parent)
{
    QScreen *screen = QGuiApplication::primaryScreen();
    QRect rec = screen->geometry();
    redraw = new QPushButton("Redraw");
    clear = new QPushButton("Clear canvas");
    setTh = new QPushButton("Set offset");
    thEdit = new QLineEdit;
    m_deltaTh = 10;

    c = new compressor;
    window = new QWidget();
    canvas = new QWidget(this);

    paintArea = new QGroupBox;
    paintArea->setTitle("Paint Area");
    paintArea->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
    controls = new QGroupBox;
    controls->setTitle("Controls");
    //controls->setSizePolicy(QSizePolicy::MinimumExpanding, QSizePolicy::MinimumExpanding);

    QHBoxLayout *h_layout_controls_1 = new QHBoxLayout;
    h_layout_controls_1->setSpacing(100);
    h_layout_controls_1->addWidget(redraw);
    h_layout_controls_1->addWidget(clear);
    QHBoxLayout *h_layout_controls_2 = new QHBoxLayout;
    h_layout_controls_2->setSpacing(100);
    h_layout_controls_2->addWidget(setTh);
    h_layout_controls_2->addWidget(thEdit);
    QVBoxLayout *v_layout_controls = new QVBoxLayout;
    v_layout_controls->addLayout(h_layout_controls_1);
    v_layout_controls->addLayout(h_layout_controls_2);
    controls->setLayout(v_layout_controls);

    QHBoxLayout *h_layout_paint = new QHBoxLayout;
    h_layout_paint->addWidget(canvas);
    paintArea->setLayout(h_layout_paint);

    QVBoxLayout *v_layout_1 = new QVBoxLayout;
    v_layout_1->addWidget(paintArea);
    v_layout_1->addWidget(controls);

    window->setLayout(v_layout_1);
    this->setGeometry(rec);
    this->setCentralWidget(window);

    //c->getRouteFromFile();

    connect(clear, SIGNAL(clicked()), this, SLOT(clearCanvas()));
    connect(redraw, SIGNAL(clicked()), this, SLOT(redrawPoly()));
    connect(setTh, SIGNAL(clicked()), this, SLOT(setOffset()));
}

void mainwindow::getArguments(int argc, char *argv[])
{
    this->argc = argc;
    this->argv = argv;
}

void mainwindow::paintEvent(QPaintEvent *)
{
    QPainter painter(this);
    QPolygonF polyF;
    QPolygonF polyFComp;
    QPen pointPen;
    QPen linePen;

    linePen.setColor(Qt::blue);
    linePen.setWidth(2);
    pointPen.setColor(Qt::blue);
    pointPen.setWidth(10);
    pointPen.setCapStyle(Qt::RoundCap);

    painter.setPen(pointPen);
    for (GroupFlight::Point &tempPoint : m_inRoute) {
        polyF << QPointF(tempPoint.x, tempPoint.y);
        painter.drawPoint(QPointF(tempPoint.x, tempPoint.y));
    }    
    painter.setPen(linePen);
    painter.drawPolygon(polyF);

    polyFComp.clear();

    pointPen.setColor(Qt::green);
    painter.setPen(pointPen);
    for (GroupFlight::Point &tempPoint1 : m_outRoute) {
        polyFComp << QPointF(tempPoint1.x, tempPoint1.y);
        painter.drawPoint(QPointF(tempPoint1.x, tempPoint1.y));
    }
    linePen.setColor(Qt::green);
    painter.setPen(linePen);
    painter.drawPolygon(polyFComp);
}

void mainwindow::mousePressEvent(QMouseEvent *event)
{
    if(event->buttons() == Qt::RightButton)
    {
        this->pointX = event->position().x();
        this->pointY = event->position().y();
        this->tempPoint.x = event->position().x();
        this->tempPoint.y = event->position().y();
        c->addValue(tempPoint);
        m_inRoute.push_back(tempPoint);
        m_outRoute = RegionReduce::parseRoute(m_inRoute, m_deltaTh);
        repaint();
    }
    else if(event->buttons() == Qt::LeftButton)
    {
        QPointF clickPoint = event->position();
        GroupFlight::Line distLine;
        distLine.p1 = GroupFlight::Point(clickPoint.x(), clickPoint.y());
        for(uint i = 0; i < m_inRoute.size(); i++)
        {
            distLine.p2 = m_inRoute.at(i);
            if(GroupFlight::lineLength(distLine) < 5)
            {
                curPointIndex = i;
                qDebug() << curPointIndex;
            }
        }
    }
}

void mainwindow::mouseMoveEvent(QMouseEvent *event)
{
    if((m_inRoute.size() >= 3) && (curPointIndex != BigNum))
    {
        QPointF clickPoint = event->position();
        m_inRoute.erase(m_inRoute.begin() + curPointIndex);
        m_inRoute.emplace(m_inRoute.begin() + curPointIndex,
                          GroupFlight::Point(clickPoint.x(), clickPoint.y()));
        m_outRoute = RegionReduce::parseRoute(m_inRoute, m_deltaTh);
        repaint();
    }
}

void mainwindow::mouseReleaseEvent(QMouseEvent *event)
{
    curPointIndex = BigNum;
}

void mainwindow::redrawPoly()
{
    m_outRoute = RegionReduce::parseRoute(m_inRoute, m_deltaTh);
}

void mainwindow::clearCanvas()
{
    m_deltaTh = 0;
    m_inRoute.clear();
    m_outRoute.clear();
    c->clearAll();
    m_deltaTh = 10;
}

void mainwindow::setOffset()
{
    c->setDeltaTh(thEdit->text().toDouble());
    m_deltaTh = thEdit->text().toDouble();
}

mainwindow::~mainwindow()
{

}
