#ifndef _SURFACE_WAVELETDECOMPOSITION_PLUGIN_H_
#define _SURFACE_WAVELETDECOMPOSITION_PLUGIN_H_

#include "plugin_interaction.h"

#include "dialog_surface_waveletDecomposition.h"
#include "decomposition.h"

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
    virtual void drawMap(View* view, MapHandlerGen* map);

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

    MapHandlerGen* initializeObject(const QString& view, const QString& filename);

private:
    Dialog_Surface_WaveletDecomposition* m_waveletDecompositionDialog;
    QAction* m_waveletDecompositionAction;

protected:
    CGoGN::Utils::ShaderColorPerVertex* m_colorPerVertexShader;
    Utils::VBO* m_positionVBO;
    Utils::VBO* m_colorVBO;

    bool m_toDraw;

    Decomposition* m_decomposition;
};

} // namespace SCHNApps

} // namespace CGoGN

#endif
