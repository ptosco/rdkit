from __future__ import annotations
import rdkit.Chem.MolDb.Loader_sa
import typing
from sqlalchemy.sql.schema import Column
from sqlalchemy.sql.sqltypes import Float
from sqlalchemy.sql.sqltypes import Integer
from sqlalchemy.sql.sqltypes import LargeBinary
from sqlalchemy.sql.sqltypes import String
from sqlalchemy.sql.sqltypes import Text
from sqlalchemy.orm.session import sessionmaker
import os
import rdkit.Chem
import rdkit.Chem.AllChem
import rdkit.Chem.Crippen
import rdkit.Chem.Descriptors
import rdkit.Chem.Lipinski
import rdkit.RDLogger
import sqlalchemy.orm.attributes
import sqlalchemy.orm.decl_api
import sqlalchemy.orm.instrumentation
import sqlalchemy.orm.mapper
import sqlalchemy.sql.schema

__all__ = [
    "AllChem",
    "Chem",
    "Column",
    "Compound",
    "ConnectToSchema",
    "Crippen",
    "Descriptors",
    "Float",
    "Integer",
    "LargeBinary",
    "Lipinski",
    "LoadDb",
    "ProcessMol",
    "RegisterSchema",
    "String",
    "Text",
    "create_engine",
    "decBase",
    "declarative_base",
    "logger",
    "logging",
    "os",
    "sessionmaker"
]


class Compound(sqlalchemy.orm.decl_api.Base):
    __mapper__: sqlalchemy.orm.mapper.Mapper # value = <Mapper at 0x7fd0588847f0; Compound>
    __table__: sqlalchemy.sql.schema.Table # value = Table('molecules', MetaData(), Column('guid', Integer(), table=<molecules>, primary_key=True, nullable=False), Column('molpkl', LargeBinary(), table=<molecules>), schema=None)
    __tablename__ = 'molecules'
    _sa_class_manager: sqlalchemy.orm.instrumentation.ClassManager # value = <ClassManager of <class 'rdkit.Chem.MolDb.Loader_sa.Compound'> at 7fd05888dd60>
    guid: sqlalchemy.orm.attributes.InstrumentedAttribute
    molpkl: sqlalchemy.orm.attributes.InstrumentedAttribute
    pass
__warningregistry__: dict # value = {'version': 74, ('Deprecated API features detected! These feature(s) are not compatible with SQLAlchemy 2.0. To prevent incompatible upgrades prior to updating applications, ensure requirements files are pinned to "sqlalchemy<2.0". Set environment variable SQLALCHEMY_WARN_20=1 to show all deprecation warnings.  Set environment variable SQLALCHEMY_SILENCE_UBER_WARNING=1 to silence this message. (Background on SQLAlchemy 2.0 at: https://sqlalche.me/e/b8d9)', <class 'sqlalchemy.exc.MovedIn20Warning'>, 20): True}
logger: rdkit.RDLogger.logger
ConnectToSchema = rdkit.Chem.MolDb.Loader_sa.RegisterSchema
decBase = sqlalchemy.orm.decl_api.Base
