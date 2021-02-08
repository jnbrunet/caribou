#include <SofaCaribou/Material/SaintVenantKirchhoffMaterial.h>
#include <SofaCaribou/Material/NeoHookeanMaterial.h>

#include <sofa/core/ObjectFactory.h>
#include <sofa/defaulttype/VecTypes.h>

using namespace sofa::core;

namespace SofaCaribou::material {

static int SaintVenantKirchhoffClass = RegisterObject("Caribou Saint-Venant-Kirchhoff hyperelastic material")
//                                                 .add< SaintVenantKirchhoffMaterial<sofa::defaulttype::Vec1Types> >()
//                                                 .add< SaintVenantKirchhoffMaterial<sofa::defaulttype::Vec2Types> >()
                                                 .add< SaintVenantKirchhoffMaterial<sofa::defaulttype::Vec3Types> >();

static int NeoHookeanClass = RegisterObject("Caribou NeoHookeanMaterial hyperelastic material")
//    .add< NeoHookeanMaterial<sofa::defaulttype::Vec1Types> >()
//    .add< NeoHookeanMaterial<sofa::defaulttype::Vec2Types> >()
    .add< NeoHookeanMaterial<sofa::defaulttype::Vec3Types> >();

} // namespace SofaCaribou::material
