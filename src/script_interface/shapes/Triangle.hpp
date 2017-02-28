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

#ifndef SCRIPT_INTERFACE_SHAPES_TRIANGLE_HPP
#define SCRIPT_INTERFACE_SHAPES_TRIANGLE_HPP

#include "Shape.hpp"
#include "core/shapes/Triangle.hpp"

namespace ScriptInterface {
namespace Shapes {


class Point : public Shape {
public:
  Point() : m_pos(new ::Shapes::Point()) {}

  const std::string name() const override { return "Shapes::Point"; }

  ParameterMap valid_parameters() const override;
  VariantMap get_parameters() const override;
  void set_parameter(const std::string &name, const Variant &value) override;

  std::shared_ptr<::Shapes::Shape> shape() const override { return m_pos; }

private:
  std::shared_ptr<::Shapes::Point> m_pos;
};


class Segment : public Shape {
public:
  Segment() : m_seg(new ::Shapes::Segment()) {}

  const std::string name() const override { return "Shapes::Segment"; }

  ParameterMap valid_parameters() const override;
  VariantMap get_parameters() const override;
  void set_parameter(const std::string &name, const Variant &value) override;

  std::shared_ptr<::Shapes::Shape> shape() const override { return m_seg; }

private:
  std::shared_ptr<::Shapes::Segment> m_seg;
};

// 
    
class Triangle : public Shape {
public:
  Triangle() : m_tri(new ::Shapes::Triangle()) {}

  const std::string name() const override { return "Shapes::Triangle"; }

  ParameterMap valid_parameters() const override;
  VariantMap get_parameters() const override;
  void set_parameter(const std::string &name, const Variant &value) override;

  std::shared_ptr<::Shapes::Shape> shape() const override { return m_tri; }

private:
  std::shared_ptr<::Shapes::Triangle> m_tri;
};

} /* namespace Shapes */
} /* namespace ScriptInterface */

#endif
