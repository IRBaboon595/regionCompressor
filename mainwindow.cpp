#include "mainwindow.h"

using namespace util;

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
    connect(this, &mainwindow::repaintSignal, [=]{this->repaint();});
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
    QPolygonF polyFCompTests;
    QPen pointPen;
    QPen linePen;

    linePen.setColor(Qt::blue);
    linePen.setWidth(2);
    pointPen.setColor(Qt::blue);
    pointPen.setWidth(10);
    pointPen.setCapStyle(Qt::RoundCap);

    painter.setPen(pointPen);
    for (GroupFlight::Point &tempPoint : m_inRoute) {
    //for (util::math::Point &tempPoint : inPoly) {
        polyF << QPointF(tempPoint.x, tempPoint.y);
        painter.drawPoint(QPointF(tempPoint.x, tempPoint.y));
    }
    painter.setPen(linePen);
    painter.drawPolygon(polyF);

    polyFComp.clear();

    pointPen.setColor(Qt::green);
    painter.setPen(pointPen);
    for (GroupFlight::Point &tempPoint1 : m_outRoute) {
    //for (util::math::Point &tempPoint1 : outPoly) {
        polyFComp << QPointF(tempPoint1.x, tempPoint1.y);
        painter.drawPoint(QPointF(tempPoint1.x, tempPoint1.y));
    }
    linePen.setColor(Qt::green);
    painter.setPen(linePen);
    painter.drawPolygon(polyFComp);

    pointPen.setColor(Qt::red);
    painter.setPen(pointPen);
    for (GroupFlight::Point &tempPoint1 : m_testRoute) {
    //for (util::math::Point &tempPoint1 : outPoly) {
        polyFCompTests << QPointF(tempPoint1.x, tempPoint1.y);
        painter.drawPoint(QPointF(tempPoint1.x, tempPoint1.y));
    }
    linePen.setColor(Qt::red);
    painter.setPen(linePen);
    painter.drawPolygon(polyFCompTests);
}

void mainwindow::mousePressEvent(QMouseEvent *event)
{
    if(event->buttons() == Qt::RightButton)
    {
        /*this->pointX = event->pos().x();
        this->pointY = event->pos().y();*/
        this->tempPoint.x = event->pos().x();
        this->tempPoint.y = event->pos().y();
        c->addValue(tempPoint);
        m_inRoute.push_back(tempPoint);
        m_outRoute = RegionReduce::parseRoute(m_inRoute, m_deltaTh, toggleVal, &m_testRoute);

        /*this->tempTestPoint.x = event->pos().x();
        this->tempTestPoint.y = event->pos().y();
        inPoly.push_back(tempTestPoint);
        outPoly = inPoly.adjusted(m_deltaTh);*/
        update();
    }
    else if(event->buttons() == Qt::LeftButton)
    {
        QPoint clickPoint = event->pos();
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

        /*util::math::Line distLine;
        distLine.p1 = util::math::Point(clickPoint.x(), clickPoint.y());

        for(uint i = 0; i < inPoly.size(); i++)
        {
            distLine.p2 = inPoly.at(i);
            if(distLine.length() < 5)
            {
                curPointIndex = i;
                qDebug() << curPointIndex;
            }
        }*/
    }
}

void mainwindow::mouseMoveEvent(QMouseEvent *event)
{
    if((m_inRoute.size() >= 3) && (curPointIndex != BigNum))
    {
        QPoint clickPoint = event->pos();
        m_inRoute.erase(m_inRoute.begin() + curPointIndex);
        m_inRoute.emplace(m_inRoute.begin() + curPointIndex,
                          GroupFlight::Point(clickPoint.x(), clickPoint.y()));

        m_outRoute = RegionReduce::parseRoute(m_inRoute, m_deltaTh, toggleVal, &m_testRoute);

        //emit repaintSignal();
    }
    update();

    /*if((inPoly.size() >= 3) && (curPointIndex != BigNum))
    {
        QPoint clickPoint = event->pos();
        inPoly.erase(inPoly.begin() + curPointIndex);
        inPoly.emplace(inPoly.begin() + curPointIndex,
                       util::math::Point(clickPoint.x(), clickPoint.y()));

        outPoly = inPoly.adjusted(m_deltaTh);

        //emit repaintSignal();
    }
    update();*/
}

void mainwindow::mouseReleaseEvent(QMouseEvent *event)
{
    curPointIndex = BigNum;
    //repaint();
}

void mainwindow::redrawPoly()
{
    //m_outRoute = RegionReduce::parseRoute(m_inRoute, m_deltaTh);
    toggleVal = toggleVal == false ? true : false;
}

void mainwindow::clearCanvas()
{
    m_deltaTh = 0;
    m_inRoute.clear();
    m_outRoute.clear();
    c->clearAll();
    m_deltaTh = 10;

    /*inPoly.clear();
    outPoly.clear();*/
}

void mainwindow::setOffset()
{
    c->setDeltaTh(thEdit->text().toDouble());
    m_deltaTh = thEdit->text().toDouble();
}

mainwindow::~mainwindow()
{

}
