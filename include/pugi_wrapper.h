#ifndef INCLUDE_PUGI_WRAPPER_H_
#define INCLUDE_PUGI_WRAPPER_H_

#include <string>
#include "errorcode.h"
#include "pugixml.hpp"

bool childExist(pugi::xml_node parent, std::string childName);

bool childAttributeExist(pugi::xml_node parent, std::string childName);

int childIntValue(pugi::xml_node parent, std::string childName, int default_value);

double childDoubleValue(pugi::xml_node parent, std::string childName, double default_value);

double childDoubleAttribute(pugi::xml_node child, std::string attribute,
		double units, bool throwIfEmpty = true);

std::string childStringValue(pugi::xml_node parent, std::string childName,
		bool throwIfEmpty = true);

std::string childStringAttribute(pugi::xml_node child, std::string attribute,
		bool throwIfEmpty = true);

pugi::xml_node childNode(pugi::xml_node parent, std::string childName,
		bool throwIfEmpty = true);

#endif /* INCLUDE_PUGI_WRAPPER_H_ */
