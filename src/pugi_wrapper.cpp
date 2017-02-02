#include "pugi_wrapper.h"

bool childExist(pugi::xml_node parent, std::string childName) {
	pugi::xml_node node = parent.child(childName.c_str());
	return (node);
}

int childIntValue(pugi::xml_node parent, std::string childName, int default_value) {
	pugi::xml_node node = parent.child(childName.c_str());
	if (!node)
		return default_value;
	return (node.attribute("value").as_int());
}

double childDoubleValue(pugi::xml_node parent, std::string childName, double default_value) {
	pugi::xml_node node = parent.child(childName.c_str());
	if (!node)
		return default_value;
	return (node.attribute("value").as_double());
}

//std::string childStringValue(xml_node parent, std::string childName, bool throwIfEmpty)
//{
//	xml_node node = parent.child(childName.c_str());
//	if (!node and throwIfEmpty) {
//		throw DRAGON::TDragonExceptions("Error reading XML card: "+childName+" not specified");
//	}
//	return (node.attribute("value").as_string());
//}
//
//bool childAttributeExist(xml_node child, std::string attribute)
//{
//	return (child.attribute(attribute.c_str()));
//}
//
//int childIntAttribute(xml_node child, std::string attribute, bool throwIfEmpty)
//{
//	if (!child and throwIfEmpty) {
//		throw DRAGON::TDragonExceptions("Error reading XML card: "+attribute+" not specified");
//	}
//	return (child.attribute(attribute.c_str()).as_int());
//}
//
//double childDoubleAttribute(xml_node child, std::string attribute, double units, bool throwIfEmpty)
//{
//	if (!child and throwIfEmpty) {
//		throw DRAGON::TDragonExceptions("Error reading XML card: "+attribute+" not specified");
//	}
//	return (child.attribute(attribute.c_str()).as_double() * units);
//}
//
//std::string childStringAttribute(xml_node child, std::string attribute, bool throwIfEmpty)
//{
//	if (!child and throwIfEmpty) {
//		throw DRAGON::TDragonExceptions("Error reading XML card: "+attribute+" not specified");
//	}
//	return (child.attribute(attribute.c_str()).as_string());
//}
//
//xml_node childNode(xml_node parent, std::string childName, bool throwIfEmpty)
//{
//	xml_node node = parent.child(childName.c_str());
//	if (!node and throwIfEmpty) {
//		throw DRAGON::TDragonExceptions("Error reading XML card: "+childName+" not specified");
//	}
//	return (node);
//}

