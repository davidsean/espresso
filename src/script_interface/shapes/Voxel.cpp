/*
  Copyright (C) 2010,2011,2012,2013,2014 The ESPResSo project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
  Max-Planck-Institute for Polymer Research, Theory Group

  This file is part of ESPResSo.

  ESPResSo is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  ESPResSo is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "Voxel.hpp"

using std::vector;
using std::string;

namespace ScriptInterface {
namespace Shapes {



ParameterMap Square::valid_parameters() const {
  return {{"a", {ParameterType::DOUBLE_VECTOR, 3, true}},
          {"b", {ParameterType::DOUBLE_VECTOR, 3, true}},
          {"c", {ParameterType::DOUBLE_VECTOR, 3, true}} };
}

VariantMap Square::get_parameters() const {
  Vector3d a;
  Vector3d b;
  int i;
  for  (i=0; i<3; i++) {
    a[i]=m_sqr->va()[i]+m_sqr->pc()[i];
    b[i]=m_sqr->vb()[i]+m_sqr->pc()[i];
  }
  return {{"a", a}, {"b", b}, {"c", m_sqr->pc()}};
}


void Square::set_parameter(const string &name,
                         const ScriptInterface::Variant &value) {
  // if C is undefined, A and B are meaningless.
  // C, or m_pc() needs to be called before m_va() and m_vb 
  // I just swapped the roles of point P and C because the corner C needs to be defined first!  
  if (name == "a") {
    /* Get the variant as vector, and explicitly construct a Vector3d
       from that. */
    m_sqr->set_pc(get_value<Vector3d>(value));
  }
  
  if (name == "b") {
    /* Get the variant as vector, and explicitly construct a Vector3d
       from that. */
    m_sqr->set_vb(get_value<Vector3d>(value));
  }
  
  if (name == "c") {
    /* Get the variant as vector, and explicitly construct a Vector3d
       from that. */
    m_sqr->set_va(get_value<Vector3d>(value));
  }
}


ParameterMap Voxel::valid_parameters() const {
  return {{"pos", {ParameterType::DOUBLE_VECTOR, 3, true}},
          {"l", {ParameterType::DOUBLE, true}}};
}

VariantMap Voxel::get_parameters() const {
  return {{"pos", m_vox->va()}, {"l", m_vox->l()}};
}


void Voxel::set_parameter(const string &name,
                         const ScriptInterface::Variant &value) {
  if (name == "pos") {
    /* Get the variant as vector, and explicitly construct a Vector3d
       from that. */
    m_vox->set_va(get_value<Vector3d>(value));
  }
  
  SET_PARAMETER_HELPER("l", m_vox->l());
    
}


} /* namespace Shapes */
} /* namespace ScriptInterface */
