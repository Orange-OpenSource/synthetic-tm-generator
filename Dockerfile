# Copyright 2018 Orange
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#  http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
FROM ubuntu:16.04
MAINTAINER Paul Chaignon <paul.chaignon@orange.com>

RUN apt-get update
RUN apt-get install -y python-pip python-tk glpk-utils
RUN pip install --upgrade pip

ADD . /workdir
WORKDIR /workdir

RUN pip install -r requirements.txt

RUN chmod a+x generate.py

ENTRYPOINT ["./generate.py"]
