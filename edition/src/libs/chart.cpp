#include "libs/chart.h"


DrawChart::DrawChart(QString title) {
  chart = new QChart();
  chart->setTitle(title);
  chart->setTitleFont(QFont("Times New Roman", 30));
}


void DrawChart::AddLine(const QVector<double>& values) {
  QVector<double> x = QVector<double>(values.size());
  std::iota(x.begin(), x.end(), 0);
  AddLine(x, values);
}

void DrawChart::AddLine(const QVector<int>& values) {
  QVector<double> y;
  y.reserve(values.size());
  std::copy(values.cbegin(), values.cend(), std::back_inserter(y));

  QVector<double> x = QVector<double>(values.size());
  std::iota(x.begin(), x.end(), 0);
  AddLine(x, y);
}

void DrawChart::AddLine(const QVector<double>& x, const QVector<double>& y) {
  switch (plot_type) {
  case 0:
    plot_type = 1;
    break;
  case 2:
    std::cerr << "Cannot mix bars and line/scatter in the same plot." << std::endl;
    return;
  default:
    break;
  }

  QLineSeries* series = new QLineSeries();
  series->setName(QString::number(series_list.size() + 1));

  int nb_values = x.size();
  if (nb_values != y.size()) {
    std::cerr << "Invalid dimension : #x != #y" << std::endl;
    return;
  }

  double min_X = x[0]; double max_X = x[0];
  double min_Y = y[0]; double max_Y = y[0];
  for (int i = 0; i < nb_values; i++) {
    series->append(x[i], y[i]);
    if (x[i] < min_X) min_X = x[i];
    if (x[i] > max_X) max_X = x[i];
    if (y[i] < min_Y) min_Y = y[i];
    if (y[i] > max_Y) max_Y = y[i];
  }
  series_list.push_back(series);

  UpdateRanges(Vector2(min_X, max_X), Vector2(min_Y, max_Y));
  chart->addSeries(series);
}


void DrawChart::AddScatter(const QVector<double>& values) {
  QVector<double> x = QVector<double>(values.size());
  std::iota(x.begin(), x.end(), 0);
  AddScatter(x, values);
}

void DrawChart::AddScatter(const QVector<int>& values) {
  QVector<double> y;
  y.reserve(values.size());
  std::copy(values.cbegin(), values.cend(), std::back_inserter(y));

  QVector<double> x = QVector<double>(values.size());
  std::iota(x.begin(), x.end(), 0);
  AddScatter(x, y);
}

void DrawChart::AddScatter(const QVector<double>& x, const QVector<double>& y) {
  switch (plot_type) {
  case 0:
    plot_type = 1;
    break;
  case 2:
    std::cerr << "Cannot mix bars and line/scatter in the same plot." << std::endl;
    return;
  default:
    break;
  }

  QScatterSeries* series = new QScatterSeries();
  series->setName(QString::number(series_list.size() + 1));

  int nb_values = x.size();
  if (nb_values != y.size()) {
    std::cerr << "Invalid dimension : #x != #y" << std::endl;
    return;
  }

  double min_X = x[0]; double max_X = x[0];
  double min_Y = y[0]; double max_Y = y[0];
  for (int i = 0; i < nb_values; i++) {
    series->append(x[i], y[i]);
    if (x[i] < min_X) min_X = x[i];
    if (x[i] > max_X) max_X = x[i];
    if (y[i] < min_Y) min_Y = y[i];
    if (y[i] > max_Y) max_Y = y[i];
  }
  series_list.push_back(series);

  UpdateRanges(Vector2(min_X, max_X), Vector2(min_Y, max_Y));
  chart->addSeries(series);
}


void DrawChart::AddBars(const QVector<int>& values) {
  switch (plot_type) {
  case 0:
    plot_type = 2;
    break;
  case 1:
    std::cerr << "Cannot mix line/scatter and bars in the same plot." << std::endl;
    return;
  default:
    break;
  }

  // Bars
  QBarSet* set0 = new QBarSet(QString::number(series_list.size() + 1));
  int nb_values = values.size();
  int max_value = 0;
  for (int i = 0; i < nb_values; i++) {
    if (values[i] > max_value) max_value = values[i];
    *set0 << values[i];
  }

  range_Y[1] = max_value;

  QBarSeries* series = new QBarSeries();
  series->append(set0);
  series_list.push_back(series);

  nb_tick_X = values.size();

  chart->addSeries(series);
}

void DrawChart::AddBars(const QVector<int>& values, double min, double max) {

  switch (plot_type) {
  case 0:
    plot_type = 2;
    break;
  case 1:
    std::cerr << "Cannot mix bars and line/scatter in the same plot." << std::endl;
    return;
  default:
    break;
  }

  // Bars
  QBarSet* set0 = new QBarSet(QString::number(series_list.size() + 1));
  nb_tick_X = values.size();
  int max_value = 0;
  for (int i = 0; i < nb_tick_X; i++) {
    if (values[i] > max_value) max_value = values[i];
    *set0 << values[i];
  }

  range_Y[1] = max_value;
  range_X = Vector2(min, max);
  float_x = true;

  QBarSeries* series = new QBarSeries();
  series->append(set0);
  series_list.push_back(series);

  chart->addSeries(series);
}


void DrawChart::UpdateRanges(Vector2 new_range_X, Vector2 new_range_Y) {
  if (series_list.size() == 0) {
    range_X = new_range_X;
    range_Y = new_range_Y;
    return;
  }

  range_X[0] = std::min(range_X[0], new_range_X[0]);
  range_X[1] = std::min(range_X[1], new_range_X[1]);
  range_Y[0] = std::min(range_Y[0], new_range_Y[0]);
  range_Y[1] = std::min(range_Y[1], new_range_Y[1]);
}


void DrawChart::SetUpChartView() {
  int nb_series = series_list.size();
  if (nb_series == 0) {
    std::cerr << "Add at least 1 series before Display()" << std::endl;
    return;
  }
  if (nb_series == 1) chart->legend()->hide();

  AddAxis();


  if (!chartView) {
    chartView = new QChartView(chart);
    chartView->setRenderHint(QPainter::Antialiasing);
    chartView->setSceneRect(QRect(QPoint(0, 0), display_window_size));
    chartView->move(250, 150);
    chartView->show();
  }
}

void DrawChart::AddAxis() {
  int nb_series = series_list.size();

  if (plot_type == 1) {

    QValueAxis* axisX = new QValueAxis;
    axisX->setRange(range_X[0], range_X[1]);
    axisX->setTickCount(nb_tick_X);
    axisX->setLabelFormat("%.2f");
    axisX->setTitleFont(QFont("Times New Roman", 30));
    axisX->setTitleBrush(QBrush(Qt::black));
    axisX->setTitleText(title_X);
    axisX->setLabelsFont(QFont("Times New Roman", 20));
    axisX->setLabelsColor(Qt::black);
    axisX->setLinePen(QPen(QBrush(QColor(0, 0, 0)), 2));
    chart->addAxis(axisX, Qt::AlignBottom);
    for (int i = 0; i < nb_series; i++)  series_list[i]->attachAxis(axisX);

    QValueAxis* axisY = new QValueAxis;
    axisY->setRange(range_Y[0], range_Y[1]);
    axisY->setTickCount(nb_tick_Y);
    axisY->setLabelFormat("%.2f");
    axisY->setTitleFont(QFont("Times New Roman", 30));
    axisY->setTitleBrush(QBrush(Qt::black));
    axisY->setTitleText(title_Y);
    axisY->setLabelsFont(QFont("Times New Roman", 20));
    axisY->setLabelsColor(Qt::black);
    axisY->setLinePen(QPen(QBrush(QColor(0, 0, 0)), 2));
    chart->addAxis(axisY, Qt::AlignLeft);
    for (int i = 0; i < nb_series; i++)  series_list[i]->attachAxis(axisY);
  }
  else {
    // X Axis
    QStringList categories;
    if (!float_x) for (int i = 0; i < nb_tick_X; i++) categories << QString::number(i);
    else for (int i = 0; i < nb_tick_X; i++) categories << QString::number(range_X[0] + (i + 0.5) * (range_X[1] - range_X[0]) / nb_tick_X, 'f', 2);
    QBarCategoryAxis* axisX = new QBarCategoryAxis();
    axisX->setTitleFont(QFont("Times New Roman", 30));
    axisX->setTitleBrush(QBrush(Qt::black));
    axisX->setTitleText(title_X);
    axisX->setLabelsFont(QFont("Times New Roman", 20));
    axisX->setLabelsColor(Qt::black);
    axisX->setLinePen(QPen(QBrush(QColor(0, 0, 0)), 2));
    axisX->append(categories);
    chart->addAxis(axisX, Qt::AlignBottom);
    for (int i = 0; i < nb_series; i++)  series_list[i]->attachAxis(axisX);
    // Y Axis
    QValueAxis* axisY = new QValueAxis();
    axisY->setRange(0, range_Y[1]);
    axisY->setLabelFormat("%i");
    axisY->setTitleFont(QFont("Times New Roman", 30));
    axisY->setTitleBrush(QBrush(Qt::black));
    axisY->setTitleText(title_Y);
    axisY->setLabelsFont(QFont("Times New Roman", 20));
    axisY->setLabelsColor(Qt::black);
    axisY->setLinePen(QPen(QBrush(QColor(0, 0, 0)), 2));
    chart->addAxis(axisY, Qt::AlignLeft);
    for (int i = 0; i < nb_series; i++)  series_list[i]->attachAxis(axisY);
  }
}

void DrawChart::Display() {
  SetUpChartView();
}

void DrawChart::ExportPDF(QString filename) {
  SetUpChartView();
  QPdfWriter writer(filename);
  writer.setResolution(160);
  writer.setPageSize(QPageSize(QSize(375, 250)));
  QPainter painter(&writer);
  chartView->render(&painter);
  painter.end();
}

void DrawChart::ExportPNG(QString filename) {
  SetUpChartView();
  QPixmap pixmap = chartView->grab();
  QFile file(filename);
  file.open(QIODevice::WriteOnly);
  pixmap.save(&file, "PNG");
}



void DrawChart::DemoChart()
{
  DrawChart graph = DrawChart("Chart Title");
  QVector<double> x_values;
  QVector<double> y_values;
  for (int i = 0; i <= 20; i++)
  {
    double x = Math::Angle(i, 20);
    x_values.push_back(x);
    y_values.push_back(sin(x));
  }
  graph.AddLine(x_values, y_values);

  graph.SetNbTickX(5);
  graph.SetNbTickX(4);

  graph.SetRangeX(0, 6);
  graph.SetRangeY(-1, 1);

  graph.SetTitleX("abscissa title");
  graph.SetTitleY("ordina title");

  graph.Display();
}