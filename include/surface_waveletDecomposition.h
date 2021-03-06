#ifndef _SURFACE_WAVELETDECOMPOSITION_PLUGIN_H_
#define _SURFACE_WAVELETDECOMPOSITION_PLUGIN_H_

#include "plugin_interaction.h"

#include "dialog_surface_waveletDecomposition.h"
#include "decomposition.h"
#include "imageCoordinates.h"

#include "camera.h"

#include "mapHandler.h"

#include "Algo/Tiling/Surface/square.h"

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

    virtual void draw(View* view){}
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

    void decomposeFromDialog();
    void saveImagesFromDialog();

public slots:

    const QString initializeObject(const QString& view, QString& filename, const bool multiple=false);
    void decompose();
    void saveImages(const QString& name);
    MapHandlerGen* drawCoarseImage(const QString& mapName);

    void project2DImageTo3DSpace(const QString& mapName);

    void triangulateMap(const QString& mapName);
    void moveUpDecomposition(const QString& mapName);
    void moveDownDecomposition(const QString& mapName);

private:
    Dialog_Surface_WaveletDecomposition* m_waveletDecompositionDialog;
    QAction* m_waveletDecompositionAction;

protected:
    Decomposition* m_decomposition;
    qglviewer::Camera* m_camera;
};

} // namespace SCHNApps

} // namespace CGoGN

#endif
