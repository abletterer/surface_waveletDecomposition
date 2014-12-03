#ifndef _DIALOG_SURFACE_WAVELETDECOMPOSITION_H_
#define _DIALOG_SURFACE_WAVELETDECOMPOSITION_H_

#include "ui_dialog_surface_waveletDecomposition.h"

namespace CGoGN
{

namespace SCHNApps
{

class SCHNApps;
class MapHandlerGen;

class Dialog_Surface_WaveletDecomposition : public QDialog, public Ui::Dialog_Surface_WaveletDecomposition
{
    Q_OBJECT

public:
    Dialog_Surface_WaveletDecomposition(SCHNApps* s);

private:
    SCHNApps* m_schnapps;
    MapHandlerGen* m_selectedMap;

public slots:
    void selectedMapChanged();
    void addMapToList(MapHandlerGen* map);
    void removeMapFromList(MapHandlerGen* map);
    void addAttributeToList(unsigned int orbit, const QString& nameAttr);
};

} // namespace SCHNApps

} // namespace CGoGN

#endif
