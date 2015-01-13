#ifndef _SURFACE_WAVELETDECOMPOSITION_PLUGIN_H_
#define _SURFACE_WAVELETDECOMPOSITION_PLUGIN_H_

#include "plugin_interaction.h"

#include "dialog_surface_waveletDecomposition.h"
#include "decomposition.h"
#include "imageCoordinates.h"

#include "camera.h"

#include "mapHandler.h"

#include "Algo/Tiling/Surface/square.h"

#include "Utils/chrono.h"

namespace CGoGN
{

namespace SCHNApps
{

class Surface_WaveletDecomposition_Plugin : public PluginInteraction
{
    Q_OBJECT
    Q_INTERFACES(CGoGN::SCHNApps::Plugin)

public:
    Surface_WaveletDecomposition_Plugin()
    {}

    ~Surface_WaveletDecomposition_Plugin()
    {}

    virtual bool enable();
    virtual void disable();

    virtual void draw(View* view);
    virtual void drawMap(View* view, MapHandlerGen* map) {}

    virtual void keyPress(View* view, QKeyEvent* event) {}
    virtual void keyRelease(View* view, QKeyEvent* event) {}
    virtual void mousePress(View* view, QMouseEvent* event) {}
    virtual void mouseRelease(View* view, QMouseEvent* event) {}
    virtual void mouseMove(View* view, QMouseEvent* event) {}
    virtual void wheelEvent(View* view, QWheelEvent* event) {}
    virtual void viewLinked(View* view) {}
    virtual void viewUnlinked(View* view) {}

private slots:
    void openWaveletDecompositionDialog();
    void closeWaveletDecompositionDialog();

    void deleteBackgroundFromDialog();

    void moveUpFromDialog();
    void moveDownFromDialog();

public slots:

    const QString initializeObject(const QString& view, QString& filename, const bool multiple = false);
    void decompose(const int max_counter = -1);
    void saveImages(const QString& name, const QString& directory = "/home/blettere/Projets/Models/Decomposition/");
    void saveAllImages(const QString& name, const QString& directory = "/home/blettere/Projets/Models/Decomposition/");
    MapHandlerGen* drawCoarseImage(const QString& mapName);

    void project2DImageTo3DSpace(const QString& mapName);
    void projectNewPointsTo3DSpace(MapHandler<PFP2>* mh_map, const std::vector<Dart>& vertices);

    void triangulateMap(const QString& mapName);
    void deleteBackground(const QString& mapName);
    void moveUpDecomposition(const QString& mapName);
    void moveDownDecomposition(const QString& mapName);

private:
    void updateDrawer();

    Dialog_Surface_WaveletDecomposition* m_waveletDecompositionDialog;
    QAction* m_waveletDecompositionAction;

protected:
    Decomposition* m_decomposition;
    qglviewer::Camera* m_camera;
    std::vector<int> m_matrix_coef;
    Utils::Drawer* m_drawer;
};

} // namespace SCHNApps

} // namespace CGoGN

#endif
