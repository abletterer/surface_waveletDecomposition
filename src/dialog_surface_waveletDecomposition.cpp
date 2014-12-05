#include "dialog_surface_waveletDecomposition.h"

#include "surface_waveletDecomposition.h"
#include "schnapps.h"
#include "mapHandler.h"

namespace CGoGN
{

namespace SCHNApps
{

Dialog_Surface_WaveletDecomposition::Dialog_Surface_WaveletDecomposition(SCHNApps* s) :
    m_schnapps(s),
    m_selectedMap(NULL)
{
    setupUi(this);

    connect(m_schnapps, SIGNAL(mapAdded(MapHandlerGen*)), this, SLOT(addMapToList(MapHandlerGen*)));
    connect(m_schnapps, SIGNAL(mapRemoved(MapHandlerGen*)), this, SLOT(removeMapFromList(MapHandlerGen*)));

    connect(list_maps, SIGNAL(itemSelectionChanged()), this, SLOT(selectedMapChanged()));

    foreach(MapHandlerGen* map, m_schnapps->getMapSet().values())
        list_maps->addItem(map->getName());
}

void Dialog_Surface_WaveletDecomposition::selectedMapChanged()
{
    QList<QListWidgetItem*> currentItems = list_maps->selectedItems();
    if(!currentItems.empty())
    {
        const QString& mapname = currentItems[0]->text();
        MapHandlerGen* mh = m_schnapps->getMap(mapname);
        m_selectedMap = mh;
    }
    else
    {
        m_selectedMap = NULL;
    }
}

void Dialog_Surface_WaveletDecomposition::addMapToList(MapHandlerGen* map)
{
    list_maps->addItem(map->getName());
}

void Dialog_Surface_WaveletDecomposition::removeMapFromList(MapHandlerGen* map)
{
    QList<QListWidgetItem*> items = list_maps->findItems(map->getName(), Qt::MatchExactly);
    if(!items.empty())
        delete items[0];

    if(m_selectedMap == map)
    {
        disconnect(m_selectedMap, SIGNAL(attributeAdded(unsigned int, const QString&)), this, SLOT(addAttributeToList(unsigned int, const QString&)));
        m_selectedMap = NULL;
    }
}

} // namespace SCHNApps

} // namespace CGoGN
